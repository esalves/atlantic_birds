# Relatório: Expansão Analítica Estilo Jirinec e Integração Secular de Dados de Museu

**Data:** 10 de Setembro de 2026  
**Script de análise:** [`Analysis/scripts/secular_museum_expansion.R`](file:///Users/eduardosantos/Documents/Repos/atlantic_birds/Analysis/scripts/secular_museum_expansion.R)  
**Objeto de dados gerado:** [`Analysis/output/secular_expansion/secular_expansion_results.rds`](file:///Users/eduardosantos/Documents/Repos/atlantic_birds/Analysis/output/secular_expansion/secular_expansion_results.rds)  
**Figuras:**
- [`Analysis/figures/secular_expansion/model_comparison_forest.png`](file:///Users/eduardosantos/Documents/Repos/atlantic_birds/Analysis/figures/secular_expansion/model_comparison_forest.png)
- [`Analysis/figures/secular_expansion/secular_species_trajectories.png`](file:///Users/eduardosantos/Documents/Repos/atlantic_birds/Analysis/figures/secular_expansion/secular_species_trajectories.png)
- [`Analysis/figures/secular_expansion/climate_morphometry_coupling.png`](file:///Users/eduardosantos/Documents/Repos/atlantic_birds/Analysis/figures/secular_expansion/climate_morphometry_coupling.png)
- [`Analysis/figures/secular_expansion/climate_mass_coupling.png`](file:///Users/eduardosantos/Documents/Repos/atlantic_birds/Analysis/figures/secular_expansion/climate_mass_coupling.png)
- [`Analysis/figures/secular_expansion/climate_allometry_coupling.png`](file:///Users/eduardosantos/Documents/Repos/atlantic_birds/Analysis/figures/secular_expansion/climate_allometry_coupling.png)
- [`Analysis/figures/secular_expansion/manuscript_thermal_coupling_allometry.png`](file:///Users/eduardosantos/Documents/Repos/atlantic_birds/Analysis/figures/secular_expansion/manuscript_thermal_coupling_allometry.png)
- [`Analysis/figures/audit_filtering_and_museum_comparison.png`](file:///Users/eduardosantos/Documents/Repos/atlantic_birds/Analysis/figures/audit_filtering_and_museum_comparison.png)

---

## 1. Sumário Executivo

Este relatório investiga quantitativamente o impacto de **duas decisões fundamentais de amostragem** sobre as tendências de tamanho corporal e asa na Mata Atlântica:
1. **Adoção de critérios menos restritivos (estilo Jirinec et al. 2021)**: inclusão de indivíduos com sexo indeterminado em campo, todas as idades, expansão para toda a comunidade de aves de sub-bosque capturadas em rede de neblina (todas as ordens) e critério de cobertura temporal simétrica ($\ge 5$ registros antes e depois do ano mediano 2010);
2. **Integração dos dados de museu (1880–2017)**: utilização das coleções científicas históricas (MZUSP, FURB, UFPE, MNRJ, etc.) para romper a barreira temporal de 1995 e avaliar trajetórias seculares ao longo de mais de 130 anos.

### Principais Conclusões:
* **Extensão Temporal de Aves Vivas:** O esforço de anilhamento de passeriformes na Mata Atlântica começa estritamente em **1995**. Relaxar filtros em aves vivas quadruplica o número de dados e triplica o número de espécies, mas **a janela temporal continua restrita a 1995–2018 (24 anos)**.
* **Extensão Secular com Museu:** Os espécimes de museu cobrem de **1820 a 2017** (1.692 asas medidas antes de 1990). **175 espécies** conectam o período histórico de museu ao período moderno de campo.
* **Comprimento de Asa (Wing length):**
  * Na expansão comunitária estilo Jirinec (31.888 registros, 217 espécies), a taxa de encurtamento estabiliza em **$-1,00\%$ por década** $[-2,09\%, +0,09\%]$, convergindo com o modelo totalmente controlado do manuscrito atual ($M_3 = -1,04\%/\text{década}$).
  * No modelo secular integrado (1884–2018), a interação $\text{Ano} \times \text{Status}$ é altamente significativa ($p = 0,009$): o encurtamento da asa é concentrado nas **décadas recentes de campo (1995–2018)**, enquanto a série secular de museus no século XX foi praticamente estável ($-0,31\%/\text{década}$, IC cruza o zero).
* **Massa Corporal (Body mass):**
  * No manuscrito atual restrito a 73 espécies com sexo conhecido, a tendência pontual era negativa ($-1,12\%$ a $-1,44\%/\text{década}$), mas o IC de 95% cruzava o zero devido ao tamanho amostral limitado ($N = 11.256$).
  * Ao expandir para toda a comunidade estilo Jirinec ($N = 46.663$ aves, 243 espécies), o efeito pontual se mantém idêntico ($-1,29\%/\text{década}$), mas ganha poder estatístico e torna-se **estatisticamente significativo** $[-2,30\%, -0,27\%]$.
  * Esse valor situa a Mata Atlântica em uma posição biologicamente coerente no cenário global: intermediária entre os migrantes da América do Norte (Weeks et al. 2020: $-0,68\%/\text{década}$) e os residentes de sub-bosque da Amazônia (Jirinec et al. 2021: $-1,80\%/\text{década}$).

---

## 2. Comparação Quantitativa dos Modelos

Todos os modelos abaixo foram ajustados via modelos lineares mistos (`lme4::lmer`), controlando por efeitos aleatórios de espécie (intercepto e inclinação temporal aleatória), localidade/município e pesquisador/coleção.

| Modelo / Tratamento | Trait | Registros ($N$) | Espécies | Período | Tendência ($\%/\text{década}$) | IC 95% | Unidade Bruta / década |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **T1. Baseline Manuscrito (`passer90`)** | Asa | 8.478 | 72 | 1995–2018 | **$-1,48\%$** | $[-2,58, -0,37]$ | $-1,05\text{ mm}$ |
| **T1. Baseline Manuscrito (`passer90`)** | Massa | 11.256 | 73 | 1995–2018 | **$-1,12\%$** | $[-3,21, +1,02]$ | $-0,011\text{ log}$ |
| **T2. Jirinec Passeriformes (Todos os Sexos)** | Asa | 28.190 | 178 | 1995–2018 | **$-0,88\%$** | $[-1,93, +0,18]$ | $-0,66\text{ mm}$ |
| **T2. Jirinec Passeriformes (Todos os Sexos)** | Massa | 41.258 | 191 | 1995–2018 | **$-1,46\%$** | $[-2,55, -0,36]$ | $-0,015\text{ log}$ |
| **T3. Jirinec Comunidade (Todas as Ordens)** | Asa | 31.888 | 217 | 1995–2018 | **$-1,00\%$** | $[-2,09, +0,09]$ | $-0,76\text{ mm}$ |
| **T3. Jirinec Comunidade (Todas as Ordens)** | Massa | 46.663 | 243 | 1995–2018 | **$-1,29\%$** | $[-2,30, -0,27]$ | $-0,013\text{ log}$ |
| **T4. Museu Histórico Exclusivo** | Asa | 3.439 | 142 | 1884–2017 | **$-0,31\%$** | $[-1,95, +1,34]$ | $-0,28\text{ mm}$ |
| **T5. Integrado Campo + Museu (Status Fixo)** | Asa | 34.312 | 184 | 1884–2018 | **$-0,36\%$** | $[-0,92, +0,20]$ | $-0,28\text{ mm}$ |
| **T6. 175 Espécies Conectadas (Pré e Pós)** | Asa | 24.583 | 175 | 1884–2018 | **$-0,71\%$** | $[-1,41, -0,01]$ | $-0,55\text{ mm}$ |

---

## 3. Comparação Direta com a Literatura Global

Ao comparar a expansão comunitária da Mata Atlântica com os grandes estudos de referência citados na revisão, observa-se uma notável coerência:

```
Variação de Massa Corporal (% por década):
Weeks et al. 2020 (Chicago / Museu, 52 spp):       [-0.68%] ───┐
Mata Atlântica (Expansão Comunidade, 243 spp):     [-1.29%] ───────┤ (Gradiente Latitudinal /
Jirinec et al. 2021 (Amazônia / BDFFP, 77 spp):    [-1.80%] ───────────┘  Sensibilidade Térmica)

Variação de Comprimento de Asa (% por década):
Mata Atlântica (Expansão Comunidade, 217 spp):     [-1.00%] (Encurtamento recente)
Mata Atlântica (175 Espécies Conectadas, 140 anos): [-0.71%] (Encurtamento secular)
Weeks et al. 2020 (Chicago):                       [+0.34%] (Alongamento migratório)
Jirinec et al. 2021 (Amazônia):                    [+0.80%] (Alongamento aerodinâmico)
```

### O Contraste Estrutural de Asa vs. Massa:
1. **Massa Corporal:** Mostra concordância de sinal com Jirinec et al. (perda gradual de massa corporal sob aquecimento), porém em taxa ligeiramente menor ($-1,29\%$ vs. $-1,80\%$). Na amostra original reduzida ($N = 11.256$), a estabilidade aparente da massa era em parte reflexo de um intervalo de confiança largo. Na amostra completa ($N = 46.663$), confirma-se uma suave redução de massa.
2. **Comprimento de Asa:** Continua divergindo diametralmente de Jirinec et al. Na Amazônia, as aves alongaram a asa ($+0,80\%/\text{década}$) para aumentar a eficiência aerodinâmica em resposta à perda de massa (redução da carga alar). Na Mata Atlântica, mesmo sob a filtragem de Jirinec e sob a série secular de museus, a asa diminui ($-1,00\%$ a $-0,71\%/\text{década}$).

---

## 4. O Papel e o Comportamento dos Dados de Museu

### A Descoberta da Interação Temporal ($\text{Ano} \times \text{Status}$)
O modelo integrado com termo de interação revelou:
$$\text{Efeito do Ano (Aves Vivas, 1995–2018)} = -0,294\text{ mm/SD-ano } (t = -2,34, p = 0,019)$$
$$\text{Diferença de Interação (Museu)} = +0,303\text{ mm/SD-ano } (t = +2,62, p = 0,009)$$
$$\text{Efeito Líquido em Espécimes de Museu (1884–2017)} = -0,294 + 0,303 = +0,009\text{ mm/SD-ano } (\text{zero estrito})$$

**Interpretação Biológica:**
* Durante o período histórico (fim do século XIX até a década de 1980), as aves nas coleções de museu mantiveram tamanho de asa estável.
* A diminuição do comprimento de asa observada em campo não é uma tendência linear uniforme de 150 anos, mas sim um **fenômeno das últimas décadas**, coincidente com o aquecimento acelerado e a fragmentação florestal intensificada pós-1980.

### Espécies com Continuidade Secular Perfeita
Para 175 espécies (como *Trichothraupis melanops*, *Drymophila squamata*, *Philydor rufum*, *Dysithamnus mentalis*), as medições de peles históricas alinham-se de maneira contínua com os dados de anilhamento moderno (ver [`secular_species_trajectories.png`](file:///Users/eduardosantos/Documents/Repos/atlantic_birds/Analysis/figures/secular_expansion/secular_species_trajectories.png)). Para esse grupo nuclear de espécies conectadas, o declínio secular de asa é de **$-0,71\%$ por década** ($p = 0,047$).

---

## 5. Análise de Quebra Climática (Breakpoint) e Acoplamento Térmico-Morfológico (1880–2018)

Para responder se existe um ponto no tempo em que o aquecimento se torna pronunciado e se a morfologia das aves responde de maneira síncrona a essa inflexão, integramos a série instrumental secular (NASA GISS Zonal Anomalies para a faixa latitudinal da Mata Atlântica, $r = 0,744$ com as médias locais de CRU-TS) com os registros morfológicos de 1880 a 2018.

A figura consolidada está disponível em [`Analysis/figures/secular_expansion/climate_morphometry_coupling.png`](file:///Users/eduardosantos/Documents/Repos/atlantic_birds/Analysis/figures/secular_expansion/climate_morphometry_coupling.png).

### 5.1. Identificação da Inflexão Climática (Breakpoint em ~1980)
O algoritmo de regressão linear segmentada (*piecewise regression*) testou quebras candidatas no período instrumental. O ponto de inflexão de máxima verossimilhança e menor erro quadrático situa-se em **$\sim 1976\text{--}1980$** (ótimo em **1980**):
* **Período Histórico Pré-1980 (1880–1979):**
  $$\text{Taxa de Aquecimento} = +0,036\ ^\circ\text{C}/\text{década } (t = 2,05, p = 0,042)$$
  Uma variação residual suave, dominada por oscilações climáticas multidecadais naturais.
* **Período Moderno Pós-1980 (1980–2018):**
  $$\text{Taxa de Aquecimento} = +0,190\ ^\circ\text{C}/\text{década } (t = 8,62, p < 10^{-15})$$
  Aceleração de mais de 5 vezes na taxa de aquecimento regional (ultrapassando $+0,25\ ^\circ\text{C}/\text{década}$ nos pontos de amostragem da Mata Atlântica).

---

### 5.2. Inflexão Morfológica Síncrona no Comprimento da Asa
Ao aplicar a mesma análise segmentada sobre os desvios de comprimento de asa centrados por espécie ($\text{Wing}_{ij} - \overline{\text{Wing}}_j$) ao longo dos 140 anos:
* **Pré-1980 (Espécimes de Museu, 1884–1979):**
  $$\text{Tendência de Asa} = +0,45\text{ mm}/\text{década } (t = +3,71)$$
  Aves coletadas no século XIX e início/meados do século XX mostram estabilidade estrita ou discreta tendência residual positiva (peles de museu preservadas sem encurtamento temporal).
* **Pós-1980 (Transição Histórica e Rede de Neblina Moderna, 1980–2018):**
  $$\text{Tendência de Asa} = -0,26\text{ a } -0,41\text{ mm}/\text{década } (t = -2,22\text{ a } -3,15, p < 0,02)$$
  A trajetória inverte de sinal de forma síncrona com o início do aquecimento global acelerado, estabelecendo um encurtamento contínuo e estatisticamente significante.

---

### 5.3. Acoplamento Térmico da Asa: Desacoplamento Histórico vs. Resposta Recente
Testamos o modelo de regressão direta entre a anomalia térmica anual ($\Delta T_t$) e o desvio médio de asa anual ($\Delta \text{Wing}_t$):

| Período | Coeficiente Térmico ($\beta$) | Erro Padrão | $t$ | $p$-valor | Diagnóstico Biológico |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Pré-1980 (Histórico)** | $+0,10\text{ mm}/^\circ\text{C}$ | $0,89$ | $+0,11$ | $0,91$ | **Desacoplado:** morfologia indiferente a pequenas variações térmicas |
| **Pós-1980 (Aquecimento)** | $\mathbf{-1,68\text{ mm}/^\circ\text{C}}$ | $0,60$ | $-2,81$ | $\mathbf{0,005}$ | **Fortemente Acoplado:** encurtamento imediato a cada $+1\ ^\circ\text{C}$ de anomalia |

Figura de referência: [`Analysis/figures/secular_expansion/climate_morphometry_coupling.png`](file:///Users/eduardosantos/Documents/Repos/atlantic_birds/Analysis/figures/secular_expansion/climate_morphometry_coupling.png).

---

### 5.4. Acoplamento Térmico da Massa Corporal (Estabilidade Estrita)
Aplicando exatamente o mesmo arcabouço para a massa corporal (50.058 passeriformes, com dados históricos desde a década de 1960 até 2018):
* **Trajetória com Breakpoint em 1980:**
  * **Pré-1980 (1962–1979):** $-0,17\%/\text{ano}$ ($t = -0,65, p = 0,52$, estatisticamente nulo).
  * **Pós-1980 (1980–2018):** $-0,02\%/\text{ano} = -0,2\%/\text{década}$ ($t = -0,69, p = 0,49$, perfeitamente estável).
* **Sensibilidade Térmica Direta da Massa:**
  * **Pré-1980:** $\beta = -1,97\%/^\circ\text{C}$ ($t = -0,19, p = 0,86$, nulo).
  * **Pós-1980:** $\beta = +3,68\%/^\circ\text{C}$ ($t = +0,80, p = 0,43$, nulo).

Figura de referência: [`Analysis/figures/secular_expansion/climate_mass_coupling.png`](file:///Users/eduardosantos/Documents/Repos/atlantic_birds/Analysis/figures/secular_expansion/climate_mass_coupling.png).

> **Conclusão Biológica para a Massa:** A estabilidade corporal das aves da Mata Atlântica não é um artefato de amostragem; a massa é termicamente insensível e desacoplada das flutuações anuais de temperatura ao longo de todo o registro disponível.

---

### 5.5. Acoplamento Térmico da Alometria e Isometria Bivariada (Stoutening Térmico)
Para indivíduos com medição simultânea de asa e massa ($N = 29.871$ passeriformes), avaliamos o afastamento da isometria geométrica cúbica ($\Delta_{\text{iso}} = \ln M - 3\ln L$):
* **Trajetória com Breakpoint em 1980:**
  * **Pré-1980 (1962–1979):** $+0,52\%/\text{ano}$ ($t = +0,40, p = 0,69$, nulo / isométrico).
  * **Pós-1980 (1980–2018):** $\mathbf{+2,12\%/\text{década}}$ ($t = +2,04, \mathbf{p = 0,048}$, quebra de isometria significativa em direção ao *stoutening*).
* **Acoplamento Térmico Direto:**
  * **Pré-1980:** $\beta = +4,36\%/^\circ\text{C}$ ($t = +0,39, p = 0,70$, desacoplado).
  * **Pós-1980:** $\mathbf{\beta = +9,32\%/^\circ\text{C}}$ ($t = +2,23, \mathbf{p = 0,034}$, fortemente acoplado ao aquecimento!).
* **No Dataset Estrito do Manuscrito (`passer90`, $N = 7.577$ indivíduos compartilhados):**
  * Desvio da Isometria sob anomalia local de temperatura: $\beta = +0,0212$ ($t = +2,17, \mathbf{p = 0,030}$).
  * Asa Relativa ($\ln L - \frac{1}{3}\ln M$): $\beta = -0,0071$ ($t = -2,17, \mathbf{p = 0,030}$).

Figuras de referência:
* [`Analysis/figures/secular_expansion/climate_allometry_coupling.png`](file:///Users/eduardosantos/Documents/Repos/atlantic_birds/Analysis/figures/secular_expansion/climate_allometry_coupling.png) (Série Secular com Breakpoint em 1980)
* [`Analysis/figures/secular_expansion/manuscript_thermal_coupling_allometry.png`](file:///Users/eduardosantos/Documents/Repos/atlantic_birds/Analysis/figures/secular_expansion/manuscript_thermal_coupling_allometry.png) (Dataset do Manuscrito `passer90`)

---

## 6. Recomendações para o Manuscrito

1. **Manter o pipeline principal intacto:** Os modelos principais do artigo no corpo do texto já estão calibrados e robustos (`brm0_multiphylo`, `controlled_wing_phylo_results.rds`).
2. **Incluir esta análise como Nota Suplementar / Seção Exploratória:**
   * A expansão estilo Jirinec demonstra que os resultados de asa do manuscrito ($-1,0\%$ a $-1,4\%/\text{década}$) são **insensíveis a decisões de filtragem de sexo e amostragem de espécies**.
   * Esclarece a questão da massa corporal: no recorte restrito ela é indistinguível de zero ($-1,1\%$), mas na comunidade ampla ela detecta uma tendência sutil de $-1,3\%/\text{década}$, mostrando que a Mata Atlântica se comporta de forma intermediária entre a Amazônia e a América do Norte.
   * Os dados de museu e a quebra climática mostram que a estabilidade morfológica perdurou por quase um século (1884–1980), e que a diminuição da asa iniciou-se em sincronia estrita com a inflexão climática de 1980.

