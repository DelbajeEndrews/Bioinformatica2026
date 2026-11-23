# Diversidade e estatística multivariada no MicrobiomeAnalyst

## O que vamos fazer

Na aula vimos como ler gráficos de composição, de diversidade alfa e de ordenação (PCoA). Agora vamos **produzir** esses gráficos e interpretá-los na plataforma online **MicrobiomeAnalyst** ([microbiomeanalyst.ca](https://www.microbiomeanalyst.ca)).

Usaremos os dados processados com o Kraken2: **comunidade do intestino do camundongo logo após o desmame (dias 0 a 9) e na fase adulta (dias 141 a 150)**

A atividade tem duas partes:

> **Atividade requer os arquivos da atividade anterior gerados pelo Kraken2?**

## Parte A — Montar as tabelas

### A.1 O que o MicrobiomeAnalyst precisa

O MicrobiomeAnalyst aceita três arquivos de texto separados por tabulação:

| Arquivo | Conteúdo | Primeira coluna |
| --- | --- | --- |
| `tabela_abundancia.txt` | número de reads de cada gênero (linhas) em cada amostra (colunas) | `#NAME` |
| `tabela_taxonomia.txt` | a linhagem de cada gênero (Kingdom, Phylum, Class, Order, Family, Genus) | `#TAXONOMY` |
| `metadados.txt` | o grupo de cada amostra | `#NAME` |

É a **tabela táxon × amostra** da aula. O banco SILVA tem vários gêneros chamados `uncultured` (em famílias diferentes), por isso cada gênero recebe um código único (`GEN_001`, `GEN_002`…), e o nome real fica na tabela de taxonomia.

A amostra Mock fica de fora: ela é um controle, não faz parte da comparação entre os grupos.

### A.2 Criar o script

Abra o terminal, ative o environment da atividade anterior e entre na pasta:

```bash
micromamba activate metagenomica
cd ~/metagenomica
ls relatorios
```

Devem aparecer os relatórios `.txt` das amostras.

Abra o **Editor de Texto**, cole o conteúdo abaixo e salve na pasta `metagenomica` com o nome **`monta_tabelas.sh`**. Não é preciso entender cada linha do `awk`: os comentários mostram o que cada bloco faz.

Obs: Diversos programas e pipelines de análises de dados metagenômicos (ex: QIIME2), por padrão, geram esses arquivos entre os resultados principais e podem ter algumas variações. Aqui usaremos um código especialmente escrito para os dados gerados pelo Kraken2.

```bash
#!/bin/bash
# Junta os relatorios do Kraken2 (so as amostras F3D*)
# em tres arquivos no formato do MicrobiomeAnalyst.

# 1) Tabela de abundancia e tabela de taxonomia
awk -F'\t' '
  FNR == 1 {
    amostra = FILENAME
    sub(/.*\//, "", amostra); sub(/\.txt$/, "", amostra)
    na++; amostras[na] = amostra
  }
  $1 ~ /g__[^|]*$/ {
    if (!($1 in id)) { ng++; id[$1] = sprintf("GEN_%03d", ng); linhagem[ng] = $1 }
    conta[id[$1], amostra] = $2
  }
  END {
    ab = "tabela_abundancia.txt"; tx = "tabela_taxonomia.txt"
    printf "#NAME" > ab
    for (j = 1; j <= na; j++) printf "\t%s", amostras[j] > ab
    printf "\n" > ab
    print "#TAXONOMY\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus" > tx
    for (i = 1; i <= ng; i++) {
      g = sprintf("GEN_%03d", i)
      printf "%s", g > ab
      for (j = 1; j <= na; j++) printf "\t%d", conta[g, amostras[j]] + 0 > ab
      printf "\n" > ab
      split("", r)
      n = split(linhagem[i], partes, "|")
      for (k = 1; k <= n; k++) r[substr(partes[k], 1, 1)] = substr(partes[k], 4)
      if (!("d" in r) && ("k" in r)) r["d"] = r["k"]
      linha = g
      split("d p c o f g", nivel, " ")
      for (k = 1; k <= 6; k++) linha = linha "\t" ((nivel[k] in r) ? r[nivel[k]] : "NA")
      print linha > tx
    }
  }' relatorios/F3D*.txt

# 2) Metadados: grupo de cada amostra
echo -e "#NAME\tGrupo" > metadados.txt
for relatorio in relatorios/F3D*.txt; do
  amostra=$(basename $relatorio .txt)
  dia=${amostra#F3D}
  if [ $dia -le 9 ]; then grupo=Dias_0a9; else grupo=Dias_141a150; fi
  echo -e "$amostra\t$grupo" >> metadados.txt
done

echo "Pronto! Arquivos criados:"
wc -l tabela_abundancia.txt tabela_taxonomia.txt metadados.txt
```

O que cada bloco faz:

- **Bloco 1:** o `awk` lê os relatórios de todas as amostras `F3D*`. Para cada linha que termina em um gênero (`g__`), guarda o número de reads daquela amostra. No fim, escreve a tabela de abundância (um gênero por linha, uma amostra por coluna, `0` quando o gênero não apareceu) e a tabela de taxonomia (a linhagem separada em colunas).
- **Bloco 2:** `${amostra#F3D}` apaga o começo `F3D` do nome e sobra só o dia (`F3D141` vira `141`). O `if` coloca as amostras com dia até 9 no grupo `Dias_0a9` e as demais no grupo `Dias_141a150`.

### A.3 Rodar e conferir

```bash
cd ~/metagenomica
bash monta_tabelas.sh
head -n 5 tabela_abundancia.txt
head -n 5 tabela_taxonomia.txt
cat metadados.txt
```
Recomendação: Abra os três arquivos .txt criados e os interprete. Observe como as informações entre eles estão conectadas.

O `metadados.txt` deve ter 19 amostras (9 do grupo `Dias_0a9` e 10 do grupo `Dias_141a150`). Copie os três arquivos para um local fácil de achar no navegador (por exemplo, a Área de Trabalho).

## Parte B — Análise no MicrobiomeAnalyst (navegador)

### B.1 Enviar os dados

1. Abra [microbiomeanalyst.ca](https://www.microbiomeanalyst.ca) e clique em **Click Here to Start**.
2. Escolha o módulo **Marker Data Profiling** (MDP), próprio para dados de gene marcador como o 16S.
3. Na página de envio, escolha o formato **texto/tabela** (*Plain text*, `.txt`) e envie:
   - **OTU/ASV Table:** `tabela_abundancia.txt`;
   - **Metadata file:** `metadados.txt`;
   - **Taxonomy table:** `tabela_taxonomia.txt`;
   - Em **Taxonomy labels:** selecione a opção `Not specific / Other`.
4. Clique em **Submit** e siga para a verificação dos dados.

### B.2 Filtragem e normalização

O site sugere filtros para remover gêneros muito raros (por exemplo, que aparecem em poucas amostras ou com poucas reads). **Mantenha os valores padrão.**

Na normalização, mantenha também o padrão e observe quais opções foram usadas. Lembre-se da aula: as amostras têm números diferentes de reads, e por isso as contagens brutas não podem ser comparadas diretamente (rarefação ou transformação em proporções deve ser aplicada).

### B.3 Composição: gráfico de barras empilhadas

1. No menu de análises, abra **Stacked Bar / Area Plot**.
2. Escolha o nível taxonômico **Phylum** e depois **Genus**, com agrupamento por amostra.
3. As figuras podem ser salvas com o botão de download da imagem no lado esquerdo da página.

Obs: Para retornar na página de seleção de resultados clique em "Analysis Overview" no topo da página.

### B.4 Diversidade alfa

1. Abra **Alpha Diversity**.
2. Escolha o índice **Shannon**, nível **Genus**, e o fator **Grupo** para separar as amostras.
3. Escolha a opção de gráfico em **box plot** e observe o resultado do teste estatístico que aparece junto (valor de p).
4. Repita com o índice **Observed** (riqueza observada).

### B.5 Diversidade beta e ordenação (PCoA)

1. Abra **Beta Diversity**.
2. Escolha o método de ordenação **PCoA**, a distância **Bray-Curtis**, o nível **Genus** e a coloração das amostras pelo fator **Grupo**.
3. Observe:
   - a porcentagem de variação explicada pelos eixos 1 e 2;
   - o resultado da **PERMANOVA** que o site mostra junto do gráfico (valor de **R²** e de **p**).
4. Repita trocando a distância para **Jaccard** (baseado em presença/ausência de táxons/feature) e compare os dois gráficos.
