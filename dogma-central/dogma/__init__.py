"""Simulador del dogma central de la biología molecular."""

from .replicacion import replicar
from .simulador import simular
from .traduccion import traducir
from .transcripcion import transcribir

__all__ = ["replicar", "transcribir", "traducir", "simular"]
