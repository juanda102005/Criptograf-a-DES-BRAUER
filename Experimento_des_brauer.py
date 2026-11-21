#!/usr/bin/env python3
"""
Experimento: semillas DES, llaves derivadas y analisis de Brauer.

Para cada una de las N semillas (por defecto 100):
  1. La semilla es una llave DES inicial de 64 bits (con paridad impar).
  2. Se derivan 16 llaves de 8 bytes con el calendario de claves de DES
     (PC-1, rotaciones de C_i y D_i, PC-1 inversa con paridad impar).
  3. Las 16 llaves son los 16 poligonos de UNA configuracion de Brauer.
  4. Se hace UN analisis de Brauer por semilla (dimension del algebra,
     dimension del centro, lazos, componentes, ...).
  5. Se construye UN grafo (quiver de Brauer) por semilla y se clasifica
     como conexo o no conexo.

Los resultados de cada semilla se guardan separados (una carpeta por semilla).

GRUPO DE CONTROL (opcional, activo por defecto): por cada semilla se genera
tambien una configuracion con 16 palabras de 8 bytes aleatorios (con paridad
impar, igual que DES). Sirve para comprobar si la estructura que aparece en
las llaves derivadas es distinta de la que aparecería por azar. Sin esta
comparacion no se puede concluir nada sobre la hipotesis.

Requiere brauer_des.py en la misma carpeta.
Dependencias: matplotlib y networkx (solo para las imagenes).
"""

import argparse
import csv
import json
import math
import os
import random
from typing import Dict, List, Tuple

from brauer_des import ConfiguracionBrauer, calendario, dibujar, hexw

# Los 128 valores de byte con cantidad impar de unos (los unicos posibles en
# una llave DES con paridad impar).
BYTES_PARIDAD_IMPAR = [b for b in range(256) if bin(b).count("1") % 2 == 1]

METRICAS = [
    "bytes_distintos",
    "bytes_valencia_1",
    "valencia_max",
    "aristas_distintas",
    "lazos",
    "componentes",
    "conexo",
    "dim_algebra",
    "dim_centro",
    "solape_consecutivo_medio",
]


# ------------------------------------------------------------------ semillas
def semilla_con_paridad(rng: random.Random) -> str:
    """Genera una llave DES de 64 bits (hex) con paridad impar en cada byte."""
    octetos = []
    for _ in range(8):
        b = rng.getrandbits(8) & 0xFE          # 7 bits de clave + bit de paridad en 0
        if bin(b).count("1") % 2 == 0:         # cantidad par de unos -> poner paridad en 1
            b |= 1
        octetos.append(b)
    return "".join(f"{b:02X}" for b in octetos)


# ------------------------------------------------------------------ analisis
def miembros_componentes(cfg: ConfiguracionBrauer) -> List[List[str]]:
    """Componentes conexas del quiver (conexidad debil), con sus poligonos."""
    padre = {w: w for w in cfg.nodos()}

    def raiz(x: str) -> str:
        while padre[x] != x:
            padre[x] = padre[padre[x]]
            x = padre[x]
        return x

    for a, b in cfg.aristas:
        padre[raiz(a)] = raiz(b)
    grupos: Dict[str, List[str]] = {}
    for w in cfg.nodos():
        grupos.setdefault(raiz(w), []).append(w)
    return sorted(grupos.values(), key=lambda g: (-len(g), g[0]))


def analizar(palabras: List[List[int]]):
    """Devuelve (configuracion, componentes, diccionario de metricas)."""
    cfg = ConfiguracionBrauer(palabras)
    comps = miembros_componentes(cfg)
    valencias = [cfg.val(v) for v in cfg.vertices]
    solapes = [len(set(a) & set(b)) for a, b in zip(palabras, palabras[1:])]
    m = {
        "bytes_distintos": len(cfg.vertices),
        "bytes_valencia_1": sum(1 for x in valencias if x == 1),
        "valencia_max": max(valencias),
        "aristas_distintas": len(cfg.aristas),
        "lazos": cfg.lazos(),
        "componentes": len(comps),
        "conexo": int(len(comps) == 1),
        "dim_algebra": cfg.dim_algebra,
        "dim_centro": cfg.dim_centro,
        "solape_consecutivo_medio": round(sum(solapes) / len(solapes), 4) if solapes else 0.0,
    }
    return cfg, comps, m


def guardar_json_semilla(ruta: str, k: int, semilla: str, palabras, cfg, comps, m):
    datos = {
        "semilla_id": k,
        "semilla_hex": semilla,
        "llaves_derivadas": {f"K{i + 1}": hexw(w) for i, w in enumerate(palabras)},
        "metricas": m,
        "clasificacion": "conexo" if m["conexo"] else "no conexo",
        "componentes": comps,
        "aristas": [f"{a} -> {b}" for a, b in
                    sorted(cfg.aristas, key=lambda e: (int(e[0][2:]), int(e[1][2:])))],
        "vertices": [
            {"byte": f"{v:02X}", "valencia": cfg.val(v), "mu": cfg.mu(v),
             "sucesores": cfg.sucesores[v]}
            for v in sorted(cfg.vertices)
        ],
    }
    with open(ruta, "w", encoding="utf-8") as f:
        json.dump(datos, f, indent=2, ensure_ascii=False)


# --------------------------------------------------------------- estadisticas
def resumen(valores: List[float]) -> Tuple[float, float, float, float]:
    n = len(valores)
    media = sum(valores) / n
    desv = math.sqrt(sum((x - media) ** 2 for x in valores) / (n - 1)) if n > 1 else 0.0
    return media, desv, min(valores), max(valores)


def welch(a: List[float], b: List[float]):
    """Estadistico t de Welch (None si no se puede calcular)."""
    if len(a) < 2 or len(b) < 2:
        return None
    ma, da, _, _ = resumen(a)
    mb, db, _, _ = resumen(b)
    se = math.sqrt(da ** 2 / len(a) + db ** 2 / len(b))
    return None if se == 0 else (ma - mb) / se


def estadisticas_globales(filas_des, filas_ctrl, ruta_txt: str, ruta_csv: str):
    lineas = []
    n = len(filas_des)
    lineas.append("=" * 100)
    lineas.append("ESTADISTICAS GLOBALES".center(100))
    lineas.append("=" * 100)
    lineas.append(f"Semillas analizadas: {n}")
    conexos = sum(r["conexo"] for r in filas_des)
    lineas.append(f"Semillas con grafo conexo (DES):     {conexos}/{n} ({100 * conexos / n:.1f}%)")
    lineas.append(f"Semillas con grafo no conexo (DES):  {n - conexos}/{n}")
    if filas_ctrl:
        cc = sum(r["conexo"] for r in filas_ctrl)
        lineas.append(f"Grupo de control con grafo conexo:   {cc}/{n} ({100 * cc / n:.1f}%)")
    lineas.append("")
    cab = f"{'metrica':<26}{'media':>9}{'desv':>9}{'min':>9}{'max':>9}"
    if filas_ctrl:
        lineas.append(f"{'':<26}{'--------- DES ---------':^36}   {'------- CONTROL -------':^36}  {'t Welch':>8}")
        lineas.append(cab + "   " + cab[26:] + f"{'':>2}")
    else:
        lineas.append(cab)

    filas_csv = []
    for met in METRICAS:
        vd = [r[met] for r in filas_des]
        md, dd, mind, maxd = resumen(vd)
        texto = f"{met:<26}{md:>9.2f}{dd:>9.2f}{mind:>9.2f}{maxd:>9.2f}"
        fila_csv = {"metrica": met, "media_DES": md, "desv_DES": dd, "min_DES": mind, "max_DES": maxd}
        if filas_ctrl:
            vc = [r[met] for r in filas_ctrl]
            mc, dc, minc, maxc = resumen(vc)
            t = welch(vd, vc)
            texto += f"   {mc:>9.2f}{dc:>9.2f}{minc:>9.2f}{maxc:>9.2f}  {('%8.2f' % t) if t is not None else '     n/a'}"
            fila_csv.update({"media_control": mc, "desv_control": dc, "min_control": minc,
                             "max_control": maxc, "t_welch": "" if t is None else t})
        lineas.append(texto)
        filas_csv.append(fila_csv)

    lineas.append("")
    lineas.append("Nota: |t| grande indica que la diferencia entre DES y el control es grande")
    lineas.append("respecto a la variabilidad entre semillas. No prueba por si sola una debilidad")
    lineas.append("criptografica; solo cuantifica la estructura del calendario de claves.")

    texto_final = "\n".join(lineas)
    print(texto_final)
    with open(ruta_txt, "w", encoding="utf-8") as f:
        f.write(texto_final + "\n")
    with open(ruta_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(filas_csv[0].keys()))
        w.writeheader()
        w.writerows(filas_csv)


# ----------------------------------------------------------------------- main
def escribir_csv(ruta: str, filas: List[dict]):
    with open(ruta, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0].keys()))
        w.writeheader()
        w.writerows(filas)


def main():
    ap = argparse.ArgumentParser(description="Experimento DES + analisis de Brauer por semilla")
    ap.add_argument("--semillas", type=int, default=100, help="numero de semillas (default 100)")
    ap.add_argument("--maestra", type=int, default=2026,
                    help="semilla del generador que produce las semillas (reproducibilidad)")
    ap.add_argument("--iteraciones", type=int, default=16, help="llaves derivadas por semilla (max 16)")
    ap.add_argument("--salida", default="resultados_experimento", help="carpeta de resultados")
    ap.add_argument("--sin-imagenes", action="store_true", help="no genera los grafos PNG")
    ap.add_argument("--sin-control", action="store_true", help="no genera el grupo de control")
    a = ap.parse_args()

    if not 1 <= a.iteraciones <= 16:
        ap.error("--iteraciones debe estar entre 1 y 16")

    os.makedirs(os.path.join(a.salida, "semillas"), exist_ok=True)

    rng = random.Random(a.maestra)
    semillas = [semilla_con_paridad(rng) for _ in range(a.semillas)]

    filas_sem, filas_llaves, filas_brauer = [], [], []
    filas_conex, filas_ctrl = [], []

    for k, semilla in enumerate(semillas, 1):
        carpeta = os.path.join(a.salida, "semillas", f"semilla_{k:03d}")
        os.makedirs(carpeta, exist_ok=True)

        # 1) llaves derivadas
        cal = calendario(semilla, a.iteraciones)
        palabras = [w for *_, w in cal]

        # 2) analisis de Brauer + conectividad
        cfg, comps, m = analizar(palabras)

        # 3) registros globales
        filas_sem.append({"semilla_id": k, "semilla_hex": semilla})
        for i, w in enumerate(palabras, 1):
            filas_llaves.append({"semilla_id": k, "semilla_hex": semilla,
                                 "iteracion": i, "K_i": hexw(w)})
        filas_brauer.append({"semilla_id": k, "semilla_hex": semilla, **m})
        filas_conex.append({
            "semilla_id": k, "semilla_hex": semilla,
            "clasificacion": "conexo" if m["conexo"] else "no conexo",
            "componentes": m["componentes"],
            "tamanos": "|".join(str(len(c)) for c in comps),
            "miembros": " | ".join("{" + ",".join(c) + "}" for c in comps),
        })

        # 4) resultados separados por semilla
        with open(os.path.join(carpeta, "llaves.txt"), "w", encoding="utf-8") as f:
            f.write(f"Semilla {k:03d}: {semilla}\n")
            for i, w in enumerate(palabras, 1):
                f.write(f"K{i} = {hexw(w)}  |  " + " ".join(f"{b:08b}" for b in w) + "\n")
        guardar_json_semilla(os.path.join(carpeta, "brauer.json"), k, semilla, palabras, cfg, comps, m)
        if not a.sin_imagenes:
            dibujar(cfg, os.path.join(carpeta, "quiver.png"))

        # 5) grupo de control: palabras aleatorias con paridad impar
        if not a.sin_control:
            rc = random.Random(f"control-{a.maestra}-{k}")
            pal_c = [[rc.choice(BYTES_PARIDAD_IMPAR) for _ in range(8)]
                     for _ in range(a.iteraciones)]
            _, _, mc = analizar(pal_c)
            filas_ctrl.append({"semilla_id": k, **mc})

        if k % 10 == 0 or k == len(semillas):
            print(f"  semilla {k}/{len(semillas)} procesada")

    escribir_csv(os.path.join(a.salida, "semillas.csv"), filas_sem)
    escribir_csv(os.path.join(a.salida, "llaves_derivadas.csv"), filas_llaves)
    escribir_csv(os.path.join(a.salida, "brauer_por_semilla.csv"), filas_brauer)
    escribir_csv(os.path.join(a.salida, "conectividad.csv"), filas_conex)
    if filas_ctrl:
        escribir_csv(os.path.join(a.salida, "brauer_control.csv"), filas_ctrl)

    print()
    estadisticas_globales(
        filas_brauer, filas_ctrl,
        os.path.join(a.salida, "estadisticas_globales.txt"),
        os.path.join(a.salida, "estadisticas_globales.csv"),
    )
    print(f"\nResultados guardados en: {a.salida}/")


if __name__ == "__main__":
    main()
