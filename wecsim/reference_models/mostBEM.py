"""Blade-element momentum loads for the published MOST IEA 15 MW turbine."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np


def _rx(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def _ry(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def _rz(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def _cross3(a, b):
    """Cross product for the three-vectors used in the MOST inner loops."""
    return np.array([
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    ])


@dataclass(frozen=True)
class MostBEM:
    """Published BEM rotor using raw IEA 15 MW AeroDyn blade tables.

    ``loads`` returns six root-force and moment components for each blade.
    The wind argument is a world-frame three-vector or a callable sampled at
    each blade element's world position.
    """

    radius: np.ndarray
    radius_interval: np.ndarray
    twist: np.ndarray
    chord: np.ndarray
    curve_angle: np.ndarray
    sweep: np.ndarray
    curve: np.ndarray
    airfoil_index: np.ndarray
    airfoil: np.ndarray
    hub_radius: float
    tip_radius: float
    rho_air: float = 1.225
    root_tol: float = 1/40
    max_iterations: int = 15
    root_func_tol: float = 1e-3
    root_tol_fallback: float = 1/40
    max_fallback_iterations: int = 50
    epsilon: float = 1e-6

    @classmethod
    def from_iea15mw(cls, blade_directory: str | Path) -> "MostBEM":
        """Build from the checked-in ``IEA-15MW`` BladeData directory."""
        directory = Path(blade_directory)
        blade = np.loadtxt(directory / "IEA-15-240-RWT_blade.dat", skiprows=6)
        radii = blade[:, 0] + 3.94  # Published hub radius, metres.
        count = int(blade[:, 6].max())
        raw_airfoils = [
            np.loadtxt(directory / f"IEA-15-240-RWT_Airfoil_{i}.dat", skiprows=54)
            for i in range(count)
        ]
        aoa_min = min(table[:, 0].min() for table in raw_airfoils)
        aoa_max = max(table[:, 0].max() for table in raw_airfoils)
        aoa = np.linspace(aoa_min, aoa_max,
                          max(len(table) for table in raw_airfoils))
        airfoil = np.empty((len(aoa), 4, count))
        for i, table in enumerate(raw_airfoils):
            angles, unique = np.unique(table[:, 0], return_index=True)
            airfoil[:, 0, i] = aoa
            for coefficient in range(1, 4):
                airfoil[:, coefficient, i] = np.interp(
                    aoa, angles, table[unique, coefficient],
                    left=np.nan, right=np.nan,
                )
        def mid(column):
            return (column[:-1] + column[1:]) / 2
        return cls(
            radius=mid(radii), radius_interval=np.diff(radii),
            twist=mid(blade[:, 4]), chord=mid(blade[:, 5]),
            curve_angle=np.deg2rad(mid(blade[:, 3])),
            sweep=mid(blade[:, 2]), curve=mid(blade[:, 1]),
            airfoil_index=blade[:-1, 6].astype(int)-1,
            airfoil=airfoil, hub_radius=float(radii[0]),
            tip_radius=float(radii[-1]),
        )

    def _coefficients(self, phi, pitch, node):
        alpha = np.rad2deg(phi-pitch) - self.twist[node]
        table = self.airfoil[:, :, self.airfoil_index[node]]
        cl, cd, cm = (
            np.interp(alpha, table[:, 0], table[:, k], left=np.nan, right=np.nan)
            for k in (1, 2, 3)
        )
        cn = cl*np.cos(phi) + cd*np.sin(phi)
        ct = cl*np.sin(phi) - cd*np.cos(phi)
        return cn, ct, cm

    def _loss(self, phi, node):
        radius = self.radius[node]
        sine = abs(np.sin(phi+1e-5))
        tip = 2/np.pi*np.arccos(np.exp(-1.5*(self.tip_radius-radius)/radius/sine))
        hub = 2/np.pi*np.arccos(np.exp(-1.5*(radius-self.hub_radius)/radius/sine))
        return tip*hub

    def _root_residual(self, phi, u_inf, pitch, omega, solidity, node,
                       negative=False):
        cn, ct, _ = self._coefficients(phi, pitch, node)
        loss = self._loss(phi, node)
        k = solidity*cn/(4*loss*np.sin(phi)**2)
        kt = solidity*ct/(4*loss*np.sin(phi)*np.cos(phi))
        if negative:
            return (np.sin(phi)*(1-k)
                    - np.cos(phi)/(omega*self.radius[node]/u_inf)*(1-kt))
        if k <= 2/3:
            a = k/(1+k)
        else:
            g1 = 2*loss*k-(10/9-loss)
            g2 = 2*loss*k-loss*(4/3-loss)
            g3 = 2*loss*k-(25/9-2*loss)+1e-6
            a = (g1-np.sqrt(g2))/g3
        at = kt/(1-kt)
        return (np.sin(phi)/(1-a)
                - np.cos(phi)/(omega*self.radius[node]/u_inf*(1+at)))

    def _bisect(self, lower, upper, function):
        left_value = function(lower)
        right_value = function(upper)
        root = upper
        iterations = 0
        while (abs((upper-lower)/lower) > self.root_tol_fallback
               and iterations < self.max_fallback_iterations
               and abs(right_value) > self.root_func_tol):
            iterations += 1
            root = (lower+upper)/2
            value = function(root)
            if left_value*value < 0:
                upper, right_value = root, value
            else:
                lower, left_value = root, value
        return root

    def _node_load(self, node, relative_wind, pitch, omega, a, at):
        velocity_ratio = relative_wind[0]/relative_wind[1]
        solidity = 3*self.chord[node]/(2*np.pi*self.radius[node])
        iteration = 1
        error = 1.0
        cn = ct = cm = np.nan
        with np.errstate(all="ignore"):
            while error > self.root_tol and iteration <= self.max_iterations:
                phi = np.arctan((1-a)/(1+at)*velocity_ratio)
                cn, ct, cm = self._coefficients(phi, pitch, node)
                loss = self._loss(phi, node)
                k = solidity*cn/(4*loss*np.sin(phi)**2)
                kt = solidity*ct/(4*loss*np.sin(phi)*np.cos(phi))
                if phi > 0 and k <= 2/3:
                    a = k/(1+k)
                elif phi > 0 and k > 2/3:
                    g1 = 2*loss*k-(10/9-loss)
                    g2 = 2*loss*k-loss*(4/3-loss)
                    g3 = 2*loss*k-(25/9-2*loss)
                    a = (g1-np.sqrt(g2))/g3
                elif phi < 0 and k > 1:
                    a = k/(k-1)
                elif k <= 1 and phi < 0:
                    a = 0
                at = kt/(1-kt)
                phi_next = np.arctan((1-a)/(1+at)*velocity_ratio)
                iteration += 1
                error = abs((phi-phi_next)/phi)

            if iteration > self.max_iterations or a == 0:
                positive = lambda phi: self._root_residual(
                    phi, relative_wind[0], pitch, omega, solidity, node)
                negative = lambda phi: self._root_residual(
                    phi, relative_wind[0], pitch, omega, solidity, node,
                    negative=True)
                if positive(self.epsilon)*positive(np.pi/2-self.epsilon) < 0:
                    phi = self._bisect(self.epsilon, np.pi/2, positive)
                elif negative(-np.pi/4)*negative(-self.epsilon) < 0:
                    phi = self._bisect(-np.pi/4, -self.epsilon, negative)
                else:
                    phi = self._bisect(np.pi/2+self.epsilon,
                                       np.pi-self.epsilon, positive)
                cn, ct, cm = self._coefficients(phi, pitch, node)
                loss = self._loss(phi, node)
                kt = solidity*ct/(4*loss*np.sin(phi)*np.cos(phi))
                at = kt/(1-kt)
                velocity_squared = (relative_wind[1]*(1+at)/np.cos(phi))**2
            else:
                velocity_squared = ((relative_wind[0]*(1-a))**2
                                    + (relative_wind[1]*(1+at))**2)
            load = (0.5*self.rho_air*self.chord[node]*velocity_squared
                    * np.array([cn, ct, cm*self.chord[node]]))
        if np.isnan([a, at]).any():
            a = at = 0
        return load, a, at

    def loads(self, hub_state, blade_pitch, wind) -> np.ndarray:
        """Return the source-convention six root loads for three blades.

        ``hub_state`` has 14 entries: world position, XYZ orientation, world
        linear and angular velocity, rotor azimuth, and rotor speed.
        """
        q = np.asarray(hub_state, dtype=float)
        if q.shape != (14,) or not np.isfinite(q).all():
            raise ValueError("hub_state must contain 14 finite values")
        if not callable(wind):
            vector = np.asarray(wind, dtype=float)
            if vector.shape != (3,):
                raise ValueError("wind must be a three-vector or callable")
            wind_at = lambda _: vector
        else:
            wind_at = wind
        hub_rotation = _rx(q[3]) @ _ry(q[4]) @ _rz(q[5]) @ _ry(np.deg2rad(6))
        rotor_velocity = q[9:12] + hub_rotation @ np.array([q[13], 0, 0])
        precone_rotation = _ry(np.deg2rad(-4))
        result = np.empty((6, 3))
        for blade in range(3):
            root_rotation = (hub_rotation @ _rx(q[12]+2*np.pi*blade/3)
                             @ precone_rotation)
            a = at = 0.0
            nodal = np.empty((len(self.radius), 3))
            for node in range(len(self.radius)):
                node_rotation = root_rotation @ _ry(self.curve_angle[node])
                position = q[:3] + root_rotation @ np.array([
                    self.curve[node], self.sweep[node], self.radius[node],
                ])
                node_velocity = q[6:9] + _cross3(rotor_velocity, position-q[:3])
                relative = node_rotation.T @ (
                    np.asarray(wind_at(position), dtype=float)-node_velocity
                )
                nodal[node], a, at = self._node_load(
                    node, relative, blade_pitch, q[13], a, at,
                )
            if np.isnan(nodal).any():
                augmented = np.vstack((np.zeros((1, 3)), nodal, np.zeros((1, 3))))
                indices = np.arange(len(augmented))
                nodal = np.column_stack([
                    np.interp(indices[1:-1], indices[np.isfinite(augmented[:, k])],
                              augmented[np.isfinite(augmented[:, k]), k])
                    for k in range(3)
                ])
            dr = self.radius_interval
            curvature = self.curve_angle
            radius_from_hub = self.radius-self.hub_radius
            result[:, blade] = [
                dr @ (nodal[:, 0]*np.cos(curvature)),
                -(dr @ nodal[:, 1]),
                dr @ (nodal[:, 0]*np.sin(-curvature)),
                dr @ (nodal[:, 1]*radius_from_hub*np.cos(np.deg2rad(4))
                      + nodal[:, 2]*np.sin(curvature)
                      + nodal[:, 0]*np.sin(-curvature)*self.sweep),
                dr @ (nodal[:, 0]*radius_from_hub),
                dr @ (nodal[:, 2]*np.cos(curvature)
                      - nodal[:, 0]*np.cos(curvature)*self.sweep
                      - nodal[:, 1]*self.curve),
            ]
        return result
