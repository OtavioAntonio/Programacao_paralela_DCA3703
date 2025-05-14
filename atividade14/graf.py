import matplotlib.pyplot as plt

# Tamanhos das mensagens em bytes
tamanhos = [8, 64, 512, 4096, 32768, 262144, 1048576]

# Resultados para 2 tarefas em 1 nó
tempo_total_1n = [0.002607, 0.002534, 0.005336, 0.010695, 0.042003, 0.177925, 1.107849]
tempo_medio_1n = [0.00000026, 0.00000025, 0.00000053, 0.00000107, 0.00000420, 0.00001779, 0.00011078]

# Resultados para 2 tarefas em 2 nós
tempo_total_2n = [0.030790, 0.036753, 0.054655, 0.118727, 0.435257, 1.975571, 7.021725]
tempo_medio_2n = [0.00000308, 0.00000368, 0.00000547, 0.00001187, 0.00004353, 0.00019756, 0.00070217]

# Configuração dos gráficos
plt.figure(figsize=(12, 5))

# Gráfico do tempo total
plt.subplot(1, 2, 1)
plt.plot(tamanhos, tempo_total_1n, label="1 nó", marker='o')
plt.plot(tamanhos, tempo_total_2n, label="2 nós", marker='s')
plt.xlabel("Tamanho da Mensagem (bytes)")
plt.ylabel("Tempo Total (s)")
plt.title("Tempo Total por Tamanho da Mensagem")
plt.xscale("log")
plt.yscale("log")
plt.grid(True)
plt.legend()

# Gráfico do tempo médio por troca
plt.subplot(1, 2, 2)
plt.plot(tamanhos, tempo_medio_1n, label="1 nó", marker='o')
plt.plot(tamanhos, tempo_medio_2n, label="2 nós", marker='s')
plt.xlabel("Tamanho da Mensagem (bytes)")
plt.ylabel("Tempo Médio por Troca (s)")
plt.title("Tempo Médio por Troca")
plt.xscale("log")
plt.yscale("log")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
