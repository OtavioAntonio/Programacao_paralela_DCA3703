import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# Dados da tabela
dados = [
    (2, 1, 1.251505, 1.240626, 1.371582),
    (2, 2, 1.213240, 1.159035, 1.221240),
    (2, 4, 0.624557, 0.580705, 0.611221),
    (2, 10, 0.302715, 0.233284, 0.247594),
    (2, 20, 0.149477, 0.119449, 0.126088),
    (4, 1, 0.714846, 0.669460, 0.723487),
    (4, 2, 0.643310, 0.596874, 0.638009),
    (4, 5, 0.266790, 0.238902, 0.260189),
    (4, 10, 0.255465, 0.118886, 0.125685),
    (8, 1, 0.348322, 0.333565, 0.384662),
    (8, 5, 5.737016, 0.126216, 0.137778),
    (8, 10, 4.698242, 2.986216, 2.695493),
]

# Rótulos do eixo X
labels = [f"{no}n-{t}t" for no, t, *_ in dados]

# Extração dos tempos
bloqueante = [b for *_, b, _, _ in dados]
nao_bloqueante = [n for *_, _, n, _ in dados]
overlap = [o for *_, _, _, o in dados]

# Plotagem
plt.figure(figsize=(12, 6))
plt.plot(labels, bloqueante, marker='o', label='Bloqueante')
plt.plot(labels, nao_bloqueante, marker='s', label='Não Bloqueante')
plt.plot(labels, overlap, marker='^', label='Overlap')

# Título e rótulos
plt.title('Comparação dos Tempos de Execução - Todas as Configurações')
plt.xlabel('Configuração (Nós-Tarefas por Nó)')
plt.ylabel('Tempo de Execução (s)')
plt.xticks(rotation=45)
plt.grid(True)
plt.legend()

# Ajuste do passo do eixo Y para 0.1s
plt.gca().yaxis.set_major_locator(MultipleLocator(0.1))

plt.tight_layout()
plt.show()
