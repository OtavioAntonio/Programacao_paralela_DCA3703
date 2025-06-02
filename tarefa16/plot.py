import matplotlib.pyplot as plt

# Dados
tamanhos = ['10x10', '20x20', '30x30', '40x40', '50x50']
num_processos = [10, 20, 30, 40, 50]
tempos = [0.010580, 0.013017, 0.020063, 0.005251, 0.008794]

plt.figure(figsize=(8,5))
plt.plot(num_processos, tempos, marker='o', linestyle='-', color='b')

plt.title('Tempo de Execução vs Número de Processos')
plt.xlabel('Número Total de Processos')
plt.ylabel('Tempo de Execução (segundos)')
plt.grid(True)

# Anotar cada ponto com o tamanho da matriz correspondente
for i, txt in enumerate(tamanhos):
    plt.annotate(txt, (num_processos[i], tempos[i]), textcoords="offset points", xytext=(0,10), ha='center')

plt.show()
