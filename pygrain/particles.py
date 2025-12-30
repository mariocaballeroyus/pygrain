from __future__ import annotations

import math
from typing import Protocol
from dataclasses import dataclass


class Particle(Protocol):

    @property
    def volume(self) -> float:
        pass


@dataclass(slots=False)
class Sphere:
    radius: float

    @property
    def volume(self) -> float:
        return (4 / 3) * math.pi * self.radius ** 3


@dataclass(slots=False)
class Spheroid:
    minor_axis: float
    aspect_ratio: float

    @property
    def volume(self) -> float:
        return (1 / 3) * math.pi * self.minor_axis**3 * self.aspect_ratio
    
    @property
    def major_axis(self) -> float:
        return self.minor_axis * self.aspect_ratio
    

@dataclass(slots=False)
class Cylinder:
    diameter: float
    aspect_ratio: float

    @property
    def volume(self) -> float:
        return (1 / 4) * math.pi * self.diameter**3 * self.aspect_ratio
    
    @property
    def length(self) -> float:
        return self.diameter * self.aspect_ratio