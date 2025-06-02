import matplotlib.pyplot as plt

# Dados
processos = ['10 (5x2)', '20 (5x4)', '30 (5x6)', '40 (5x8)', '50 (5x10)']
tempos = [0.028930, 0.023207, 0.018615, 0.023459, 0.039994]
matrizes = ['10x10', '20x20', '30x30', '40x40', '50x50']

# Plot
plt.figure(figsize=(10, 6))
plt.plot(processos, tempos, marker='o', linestyle='-', color='blue', label='Tempo de Execução')
for i, txt in enumerate(matrizes):
    plt.annotate(txt, (processos[i], tempos[i] + 0.001), ha='center', fontsize=9)

# Estética
plt.title('Tempo de Execução vs Número de Processos MPI')
plt.xlabel('Total de Processos (Nodes x Tarefas por Node)')
plt.ylabel('Tempo de Execução (segundos)')
plt.grid(True)
plt.legend()
plt.tight_layout()

# Exibir
plt.show()
