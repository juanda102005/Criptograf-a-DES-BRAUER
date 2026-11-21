#!/usr/bin/env python3
"""
Analisis de Brauer sobre el calendario de claves DES.

Para cada clave DES de 64 bits:
  1. Se calcula el calendario de claves (PC-1, rotaciones de C_i y D_i).
  2. En cada iteracion i se reconstruye una palabra K_i de 8 bytes
     (PC-1 inversa + bit de paridad impar en cada byte).
  3. Las palabras K_1..K_i son los poligonos w_1..w_i de una configuracion de
     Brauer; los vertices son los valores de byte (00..FF).
  4. Se calculan: dimension del algebra, dimension del centro y
     conexidad del quiver (grafo de cubrimiento).

Convenciones (basadas en Partitura.py):
  - Orden ciclico de poligonos alrededor de un vertice = orden de iteracion;
    un byte repetido dentro de una misma palabra cuenta tantas veces como aparece.
  - Flechas: pares consecutivos de la secuencia ciclica de cada vertice.
  - m(v) = 2 si val(v) = 1, y m(v) = 1 en los demas casos.
  - dim = 2n + sum_v val(v) * (m(v)*val(v) - 1),  n = numero de poligonos.
  - dim Z = 1 + n + 2*(poligonos con lazo)   (formula de Partitura.py).

Si VALENCIA_COMO_PARTITURA = True se reproduce el comportamiento original de
Partitura.py, donde la valencia incluye el elemento de cierre (apariciones + 1).
"""

import argparse
import csv
import os
import random
from typing import Dict, List, Set, Tuple

VALENCIA_COMO_PARTITURA = False

# ---------------------------------------------------------------- tablas DES
PC1_C = [57, 49, 41, 33, 25, 17, 9,
         1, 58, 50, 42, 34, 26, 18,
         10, 2, 59, 51, 43, 35, 27,
         19, 11, 3, 60, 52, 44, 36]
PC1_D = [63, 55, 47, 39, 31, 23, 15,
         7, 62, 54, 46, 38, 30, 22,
         14, 6, 61, 53, 45, 37, 29,
         21, 13, 5, 28, 20, 12, 4]
PC1 = PC1_C + PC1_D                      # 56 posiciones de la clave de 64 bits
ROTACIONES = [1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1]


# ------------------------------------------------------------- utilidades bits
def hex_a_bits(h: str) -> List[int]:
    h = h.replace(" ", "")
    if len(h) != 16:
        raise ValueError("La clave debe tener 16 simbolos hexadecimales (64 bits).")
    return [int(b) for b in bin(int(h, 16))[2:].zfill(64)]


def bits_a_bytes(bits: List[int]) -> List[int]:
    return [int("".join(map(str, bits[i:i + 8])), 2) for i in range(0, 64, 8)]


def rotar_izq(bits: List[int], k: int) -> List[int]:
    return bits[k:] + bits[:k]


def bits_a_str(bits: List[int], grupo: int = 7) -> str:
    s = "".join(map(str, bits))
    return " ".join(s[i:i + grupo] for i in range(0, len(s), grupo))


# ------------------------------------------------------- calendario de claves
def pc1(clave64: List[int]) -> Tuple[List[int], List[int]]:
    """Aplica PC-1: devuelve (C0, D0). Los bits 8,16,...,64 se descartan."""
    return [clave64[p - 1] for p in PC1_C], [clave64[p - 1] for p in PC1_D]


def pc1_inversa(c: List[int], d: List[int]) -> List[int]:
    """Reconstruye 64 bits desde C_i D_i y fija paridad impar en cada byte."""
    registro = c + d
    word = [0] * 64
    for j, pos in enumerate(PC1):
        word[pos - 1] = registro[j]
    for b in range(8):                     # bit de paridad = posicion 8, 16, ...
        unos = sum(word[8 * b: 8 * b + 7])
        word[8 * b + 7] = 0 if unos % 2 == 1 else 1
    return word


def calendario(clave_hex: str, iteraciones: int = 16):
    """Devuelve lista de (i, C_i, D_i, K_i como 8 bytes)."""
    c, d = pc1(hex_a_bits(clave_hex))
    salida = []
    for i in range(1, iteraciones + 1):
        c = rotar_izq(c, ROTACIONES[i - 1])
        d = rotar_izq(d, ROTACIONES[i - 1])
        salida.append((i, c, d, bits_a_bytes(pc1_inversa(c, d))))
    return salida


# ------------------------------------------------------------ analisis Brauer
class ConfiguracionBrauer:
    """Poligonos = palabras (listas de bytes). Replica la logica de Partitura."""

    def __init__(self, palabras: List[List[int]]):
        self.datos = palabras
        self.n = len(palabras)
        self.vertices: Set[int] = {b for w in palabras for b in w}
        self.sucesores: Dict[int, List[str]] = {}
        self.aristas: Set[Tuple[str, str]] = set()
        self.dim_algebra = 0
        self.dim_centro = 0
        self.analizar()
    def lazos_flechas(self) -> int:
        """Flechas del quiver que son lazos (w_i -> w_i), una por cada par
        consecutivo igual en la secuencia cíclica de sucesores de cada byte."""
        return sum(
            1
            for seq in self.sucesores.values()
            for a, b in zip(seq, seq[1:])
            if a == b
        )
    def analizar(self):
        for v in self.vertices:
            seq = []
            for i, w in enumerate(self.datos):
                seq.extend([f"w_{i + 1}"] * w.count(v))
            seq.append(seq[0])             # cierre del ciclo
            self.sucesores[v] = seq
        for seq in self.sucesores.values():
            for a, b in zip(seq, seq[1:]):
                self.aristas.add((a, b))
        self.dim_algebra = sum(
            self.val(v) * (self.mu(v) * self.val(v) - 1) for v in self.vertices
        ) + 2 * self.n
        n_C = sum(1 for v in self.vertices if self.val(v) == 1 and self.mu(v) > 1)
        self.dim_centro = (1 + sum(self.mu(v) for v in self.vertices) + self.n
                           - len(self.vertices) + self.lazos_flechas() - n_C)

    def val(self, v: int) -> int:
        largo = len(self.sucesores[v])
        return largo if VALENCIA_COMO_PARTITURA else largo - 1

    def mu(self, v: int) -> int:
        return 2 if self.val(v) <= 1 else 1

    def nodos(self) -> List[str]:
        return [f"w_{i + 1}" for i in range(self.n)]

    def lazos(self) -> int:
        return sum(1 for w in self.nodos() if (w, w) in self.aristas)

    def componentes(self) -> int:
        padre = {w: w for w in self.nodos()}

        def raiz(x):
            while padre[x] != x:
                padre[x] = padre[padre[x]]
                x = padre[x]
            return x

        for a, b in self.aristas:
            padre[raiz(a)] = raiz(b)
        return len({raiz(w) for w in self.nodos()})

    def conexo(self) -> bool:
        return self.componentes() == 1


# ------------------------------------------------------------------ reportes
def hexw(w: List[int]) -> str:
    return " ".join(f"{b:02X}" for b in w)


def reporte_detalle(clave_hex: str, iteraciones: int):
    """Imprime todo el procedimiento para una clave (para comparar con papel y lapiz)."""
    bits = hex_a_bits(clave_hex)
    print("CLAVE:", " ".join(clave_hex[i:i + 2] for i in range(0, 16, 2)).upper())
    print("BINARIO:", " ".join("".join(map(str, bits[i:i + 8])) for i in range(0, 64, 8)))
    c0, d0 = pc1(bits)
    print("C0 =", bits_a_str(c0))
    print("D0 =", bits_a_str(d0))
    cal = calendario(clave_hex, iteraciones)
    for i, c, d, w in cal:
        print(f"\nIteracion {i} (rotacion {ROTACIONES[i - 1]})")
        print("  C%d =" % i, bits_a_str(c))
        print("  D%d =" % i, bits_a_str(d))
        print("  K%d =" % i, hexw(w), "|", " ".join(f"{b:08b}" for b in w))
    print("\nANALISIS DE BRAUER POR ITERACION")
    print(f"{'i':>3} {'vert':>5} {'dim':>6} {'dimZ':>6} {'lazos':>6} {'comp':>5} conexo")
    for i in range(1, iteraciones + 1):
        cfg = ConfiguracionBrauer([w for *_, w in cal[:i]])
        print(f"{i:>3} {len(cfg.vertices):>5} {cfg.dim_algebra:>6} {cfg.dim_centro:>6} "
              f"{cfg.lazos():>6} {cfg.componentes():>5} {'si' if cfg.conexo() else 'no'}")
    print("\nSECUENCIAS DE SUCESORES (iteracion final)")
    cfg = ConfiguracionBrauer([w for *_, w in cal])
    for v in sorted(cfg.vertices):
        print(f"  {v:02X} | val={cfg.val(v)} m={cfg.mu(v)} | {' < '.join(cfg.sucesores[v])}")


def dibujar(cfg: ConfiguracionBrauer, ruta: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import networkx as nx
    G = nx.DiGraph()
    G.add_nodes_from(cfg.nodos())
    G.add_edges_from(cfg.aristas)
    pos = nx.circular_layout(G)
    plt.figure(figsize=(8, 6))
    nx.draw_networkx_nodes(G, pos, node_color="white", edgecolors="black",
                           node_size=900, linewidths=2)
    nx.draw_networkx_edges(G, pos, edge_color="black", width=1.5, arrowsize=15,
                           connectionstyle="arc3,rad=0.1", node_size=900)
    nx.draw_networkx_labels(G, pos, labels={n: f"${n}$" for n in G.nodes()}, font_size=12)
    plt.axis("off")
    plt.margins(0.15)
    plt.savefig(ruta, dpi=300, bbox_inches="tight")
    plt.close()


def generar_claves(n: int, semilla: int) -> List[str]:
    rng = random.Random(semilla)
    return [f"{rng.getrandbits(64):016X}" for _ in range(n)]


def main():
    ap = argparse.ArgumentParser(description="Analisis de Brauer del calendario de claves DES")
    ap.add_argument("--claves", type=int, default=100, help="numero de claves aleatorias")
    ap.add_argument("--semilla", type=int, default=2026)
    ap.add_argument("--iteraciones", type=int, default=16)
    ap.add_argument("--clave", help="una clave concreta en hex (16 simbolos)")
    ap.add_argument("--detalle", action="store_true", help="reporte paso a paso (usar con --clave)")
    ap.add_argument("--dibujar", action="store_true", help="dibuja el quiver de la ultima iteracion")
    ap.add_argument("--imagenes", metavar="CARPETA",
                    help="guarda un quiver (iteracion final) por clave en esta carpeta")
    ap.add_argument("--salida", default="resultados_brauer_des.csv")
    a = ap.parse_args()

    claves = [a.clave.replace(" ", "").upper()] if a.clave else generar_claves(a.claves, a.semilla)

    if a.detalle:
        reporte_detalle(claves[0], a.iteraciones)
        if a.dibujar:
            cal = calendario(claves[0], a.iteraciones)
            dibujar(ConfiguracionBrauer([w for *_, w in cal]), "Grafo_Brauer_DES.png")
        return

    filas = []
    if a.imagenes:
        os.makedirs(a.imagenes, exist_ok=True)
    for k, clave in enumerate(claves, 1):
        cal = calendario(clave, a.iteraciones)
        if a.imagenes:
            cfg_final = ConfiguracionBrauer([w for *_, w in cal])
            dibujar(cfg_final, os.path.join(a.imagenes, f"Quiver_clave_{k:03d}_{clave}.png"))
        for i in range(1, a.iteraciones + 1):
            cfg = ConfiguracionBrauer([w for *_, w in cal[:i]])
            filas.append({
                "clave_n": k, "clave_hex": clave, "iteracion": i,
                "K_i": hexw(cal[i - 1][3]),
                "vertices": len(cfg.vertices), "dim_algebra": cfg.dim_algebra,
                "dim_centro": cfg.dim_centro, "lazos": cfg.lazos(),
                "componentes": cfg.componentes(), "conexo": int(cfg.conexo()),
            })
    with open(a.salida, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0].keys()))
        w.writeheader()
        w.writerows(filas)

    print(f"{len(claves)} claves x {a.iteraciones} iteraciones -> {a.salida}")
    print(f"{'i':>3} {'dim media':>10} {'dimZ media':>11} {'% conexo':>9} {'comp media':>11}")
    for i in range(1, a.iteraciones + 1):
        fi = [r for r in filas if r["iteracion"] == i]
        print(f"{i:>3} {sum(r['dim_algebra'] for r in fi) / len(fi):>10.1f} "
              f"{sum(r['dim_centro'] for r in fi) / len(fi):>11.1f} "
              f"{100 * sum(r['conexo'] for r in fi) / len(fi):>8.1f}% "
              f"{sum(r['componentes'] for r in fi) / len(fi):>11.2f}")


if __name__ == "__main__":
    main()

