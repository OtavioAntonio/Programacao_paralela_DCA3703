import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import re

# Função para extrair o tempo do nome do arquivo
def extrair_tempo(nome):
    match = re.search(r't(\d+)', nome)
    return int(match.group(1)) if match else -1

# Ordena os arquivos por tempo
arquivos = sorted(glob.glob("velocidade_u_fatia_z_t*.csv"), key=extrair_tempo)

# Seleciona até 6 arquivos uniformemente espaçados
num_graficos = 6
arquivos_amostrados = [arquivos[i] for i in np.linspace(0, len(arquivos) - 1, num_graficos, dtype=int)]

cols = 3
rows = 2

# Figura e eixos
fig, axs = plt.subplots(rows, cols, figsize=(14, 8), constrained_layout=True)
axs = np.array(axs).reshape(rows, cols)

# Para a barra de cores
im = None

for idx, arq in enumerate(arquivos_amostrados):
    dados = np.loadtxt(arq, delimiter=',')
    r, c = divmod(idx, cols)
    ax = axs[r][c]
    im = ax.imshow(dados, cmap='viridis', origin='lower')
    ax.set_title(os.path.basename(arq), fontsize=10)
    ax.axis('off')

# Remove subplots vazios
for i in range(len(arquivos_amostrados), rows * cols):
    r, c = divmod(i, cols)
    axs[r][c].axis('off')

# Adiciona uma barra de cores ao lado direito da figura
cbar = fig.colorbar(im, ax=axs, location='right', shrink=0.6, label='Velocidade')

plt.show()
