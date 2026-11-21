# Análisis de Brauer sobre el calendario de claves DES
 
Experimento de criptografía: ¿qué estructura aparece cuando se estudian con **configuraciones de Brauer** las 16 subllaves que DES deriva de una misma llave inicial?
 
## Qué se hizo
 
1. Se generaron **100 semillas**: llaves DES de 64 bits con paridad impar.
2. De cada semilla, el calendario de claves de DES (PC-1, rotaciones de `C_i` y `D_i`, PC-1 inversa con paridad impar) produce **16 subllaves** de 8 bytes. En total, 1600 subllaves.
3. Las 16 subllaves de una semilla forman **una** configuración de Brauer:
   - cada subllave `K_i` es un **polígono** (nodo `w_i` del quiver),
   - cada **valor de byte** (00–FF) es un **vértice**,
   - el orden cíclico alrededor de un vértice es el orden de las iteraciones, y los bytes repetidos dentro de una subllave cuentan cada vez que aparecen.
4. Por semilla se calcula: dimensión del álgebra, dimensión del centro, lazos, componentes y conexidad del quiver.
5. Se añadió un **grupo de control**: por cada semilla, 16 palabras de 8 bytes aleatorios (con paridad impar), analizadas igual, para saber si lo observado difiere del azar.
Convenciones: valencia = número de apariciones del byte; multiplicidad `m = 2` si la valencia es 1 y `m = 1` en los demás casos; `dim = 2n + Σ val(v)·(m(v)·val(v) − 1)`, con `n` el número de polígonos.
 
Dimensión del centro: se usa la fórmula de Sierra (ver Referencias),
 
`dim Z = 1 + Σ m(v) + n − |Γ₀| + #lazos − |C|`,
 
donde `#lazos` es el número de flechas del quiver que son lazos (`w_i → w_i`) y `C` es el conjunto de vértices con valencia 1 y `m > 1`. Con las convenciones de arriba se simplifica a `dim Z = 1 + n + #lazos`.
 
## Archivos
 
| Archivo | Qué es |
|---|---|
| `brauer_des.py` | Módulo base: calendario de claves DES, construcción de la configuración de Brauer, cálculos y dibujo del quiver. También se puede ejecutar solo (`--detalle` muestra el procedimiento paso a paso para una clave). |
| `experimento_des_brauer.py` | Script del experimento completo (semillas, análisis, grafos, conectividad, control y estadísticas). Importa `brauer_des.py`. |
| `resultados_experimento/` | Carpeta con todos los resultados (ver abajo). |
 
## Resultados (`resultados_experimento/`)
 
| Archivo / carpeta | Contenido |
|---|---|
| `semillas.csv` | Las 100 semillas (`semilla_id`, `semilla_hex`). |
| `llaves_derivadas.csv` | Las 1600 subllaves (16 por semilla), en hexadecimal. |
| `brauer_por_semilla.csv` | Una fila por semilla con las métricas del análisis de Brauer. |
| `brauer_control.csv` | Las mismas métricas para el grupo de control. |
| `conectividad.csv` | Clasificación de cada grafo (conexo / no conexo), número de componentes y sus miembros. |
| `estadisticas_globales.txt` / `.csv` | Media, desviación, mínimo, máximo y t de Welch, DES frente a control. |
| **`semillas/semilla_001/` … `semilla_100/`** | **Una carpeta por semilla (resultados separados).** |
 
### Los 100 grafos
 
Están en `resultados_experimento/semillas/semilla_XXX/quiver.png`, uno por semilla (de `semilla_001` a `semilla_100`). Cada imagen es el **quiver de Brauer** de esa semilla: 16 nodos (`w_1` … `w_16`, es decir, las subllaves `K_1` … `K_16`) y una flecha entre subllaves consecutivas que comparten un byte. Un lazo (`w_i → w_i`) aparece cuando un byte está en una sola subllave o repetido dentro de ella.
 
Cada carpeta `semilla_XXX/` contiene además:
 
- `llaves.txt`: las 16 subllaves en hexadecimal y binario.
- `brauer.json`: métricas, vértices (con valencia y sucesores), aristas y componentes.
- `quiver.png`: el grafo.
## Resultados principales
 
Los 100 grafos de DES resultaron **conexos**, y también los 100 del control, por lo que la conexidad no distingue entre ambos. Lo que sí los separa (medias sobre las 100 semillas):
 
| Métrica | DES | Control |
|---|---|---|
| Bytes distintos (de 128 posibles) | 35.3 | 80.9 |
| Bytes con valencia 1 | 0.0 | 46.9 |
| Polígonos con al menos un lazo (de 16) | 3.3 | 15.7 |
| Dimensión del álgebra | 496.1 | 206.7 |
| Dimensión del centro | 20.5 | 69.0 |
| Solape medio entre subllaves consecutivas (bytes) | 5.29 | 0.47 |
 
En DES ningún byte tiene valencia 1 (todos aparecen en al menos dos subllaves), mientras que en el control casi la mitad de los bytes aparecen una sola vez y cada uno produce un lazo. Esto explica que la dimensión del centro sea mucho mayor en el control y la del álgebra mayor en DES (valencias más altas).
 
Las subllaves de DES salen de rotar los mismos registros, así que comparten muchos bytes por construcción. El experimento **cuantifica esa estructura**; no demuestra por sí solo una debilidad criptográfica de DES.
 
## Cómo reproducirlo
 
```bash
python -m pip install matplotlib networkx
python experimento_des_brauer.py
```
 
Los archivos `brauer_des.py` y `experimento_des_brauer.py` deben estar en la misma carpeta. Opciones útiles:
 
- `--semillas N`: número de semillas (100 por defecto).
- `--maestra 2026`: semilla del generador; con el mismo valor el experimento es reproducible.
- `--salida carpeta`: carpeta de resultados (por defecto `resultados_experimento`).
- `--sin-imagenes`: omite los PNG.
- `--sin-control`: omite el grupo de control.
## Referencias
 
- A. Sierra, *The Dimension of the Center of a Brauer Configuration Algebra*, arXiv:1707.08562.
- E. L. Green y S. Schroll, *Brauer Configuration Algebras: A Generalization of Brauer Graph Algebras*, Bull. Sci. Math. 141 (2017) 539–572.
 