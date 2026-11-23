# Diversidade e estatística multivariada no MicrobiomeAnalyst: roteiro passo a passo

## O que vamos fazer

Na aula vimos como ler gráficos de composição, de diversidade alfa e de ordenação (PCoA). Agora vamos **produzir** esses gráficos e interpretá-los com uma ferramenta online, o **MicrobiomeAnalyst** ([microbiomeanalyst.ca](https://www.microbiomeanalyst.ca)), sem precisar programar a parte estatística.

A pergunta é a mesma da atividade com o Kraken2: **a comunidade do intestino do camundongo é diferente logo após o desmame (dias 0 a 9) e na fase adulta (dias 141 a 150)?**

A atividade tem duas partes:

| Parte | Onde | O que faz |
| --- | --- | --- |
| A | terminal | junta os relatórios do Kraken2 em três tabelas no formato do MicrobiomeAnalyst |
| B | navegador | envia as tabelas e faz composição, diversidade alfa, PCoA e PERMANOVA |

> **Não fez a atividade com o Kraken2?** Pule a parte A. Na parte B, use um dos **conjuntos de dados de exemplo** que o próprio MicrobiomeAnalyst oferece na página de envio de dados, e responda às mesmas perguntas para esse conjunto.

## Parte A — Montar as tabelas (terminal)

### A.1 O que o MicrobiomeAnalyst precisa

O MicrobiomeAnalyst aceita três arquivos de texto separados por tabulação:

| Arquivo | Conteúdo | Primeira coluna |
| --- | --- | --- |
| `tabela_abundancia.txt` | número de reads de cada gênero (linhas) em cada amostra (colunas) | `#NAME` |
| `tabela_taxonomia.txt` | a linhagem de cada gênero (Kingdom, Phylum, Class, Order, Family, Genus) | `#TAXONOMY` |
| `metadados.txt` | o grupo de cada amostra | `#NAME` |

É a **tabela táxon × amostra** da aula. Como o banco SILVA tem vários gêneros chamados `uncultured` (em famílias diferentes), cada gênero recebe um código único (`GEN_001`, `GEN_002`…), e o nome real fica na tabela de taxonomia.

A amostra Mock fica de fora: ela é um controle, não faz parte da comparação entre os grupos.

### A.2 Criar o script

Abra o terminal, ative o ambiente da atividade anterior e entre na pasta:

```bash
micromamba activate metagenomica
cd ~/metagenomica
ls relatorios
```

Devem aparecer os relatórios `.txt` das 20 amostras.

Abra o **Editor de Texto**, cole o conteúdo abaixo e salve na pasta `metagenomica` com o nome **`monta_tabelas.sh`**. Não é preciso entender cada linha do `awk`: os comentários mostram o que cada bloco faz.

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

O `metadados.txt` deve ter 19 amostras (9 do grupo `Dias_0a9` e 10 do grupo `Dias_141a150`). Copie os três arquivos para um local fácil de achar no navegador (por exemplo, a Área de Trabalho).

## Parte B — Análise no MicrobiomeAnalyst (navegador)

### B.1 Enviar os dados

1. Abra [microbiomeanalyst.ca](https://www.microbiomeanalyst.ca) e clique em **Start Here** (ou **Get Started**).
2. Escolha o módulo **Marker Data Profiling** (MDP), próprio para dados de gene marcador como o 16S.
3. Na página de envio, escolha o formato **texto/tabela** (*Plain text*, `.txt`) e envie:
   - **Taxonomic profile / abundance:** `tabela_abundancia.txt`;
   - **Metadata file:** `metadados.txt`;
   - **Taxonomy table:** `tabela_taxonomia.txt` (indique que a taxonomia está em um **arquivo separado**, e não dentro da tabela de abundância).
4. Clique em **Submit** e siga para a verificação dos dados. Confira se aparecem **19 amostras** e o fator **Grupo** com dois níveis.

### B.2 Filtragem e normalização

O site sugere filtros para remover gêneros muito raros (por exemplo, que aparecem em poucas amostras ou com poucas reads). **Mantenha os valores padrão.**

Na normalização, mantenha também o padrão e anote quais opções foram usadas. Lembre-se da aula: as amostras têm números diferentes de reads, e por isso as contagens brutas não podem ser comparadas diretamente (rarefação ou transformação em proporções).

### B.3 Composição: gráfico de barras empilhadas

1. No menu de análises, abra **Stacked Bar / Area Plot**.
2. Escolha o nível taxonômico **Phylum** e depois **Genus**, com agrupamento por amostra.
3. Salve as duas figuras (botão de download da imagem).

### B.4 Diversidade alfa

1. Abra **Alpha Diversity**.
2. Escolha o índice **Shannon**, nível **Genus**, e o fator **Grupo** para separar as amostras.
3. Escolha a opção de gráfico em **box plot** e anote o resultado do teste estatístico que aparece junto (valor de p).
4. Repita com o índice **Observed** (riqueza observada).

### B.5 Diversidade beta e ordenação (PCoA)

1. Abra **Beta Diversity**.
2. Escolha o método de ordenação **PCoA**, a distância **Bray-Curtis**, o nível **Genus** e colora as amostras pelo fator **Grupo**.
3. Anote:
   - a porcentagem de variação explicada pelos eixos 1 e 2;
   - o resultado da **PERMANOVA** que o site mostra junto do gráfico (valor de **R²** e de **p**).
4. Repita trocando a distância para **Jaccard** (presença/ausência) e compare os dois gráficos.

### B.6 Comparar com o terminal (opcional)

Se você fez a atividade com o Kraken2, compare o Shannon calculado pelo site com o do arquivo `diversidade_alfa.txt`. Os valores não precisam ser idênticos: o site filtra gêneros raros e normaliza os dados antes do cálculo.

## Perguntas para o relatório

1. Quais filos dominam as amostras? Algum filo muda claramente de proporção entre os dois grupos?
2. Qual grupo tem maior diversidade alfa (Shannon e riqueza observada)? A diferença foi estatisticamente significativa?
3. No gráfico de PCoA (Bray-Curtis), as amostras dos dois grupos se separam? Quanto da variação os dois primeiros eixos explicam juntos?
4. Interprete o R² e o p da PERMANOVA com suas próprias palavras.
5. O padrão da PCoA mudou ao usar a distância de Jaccard em vez de Bray-Curtis? O que essa diferença diz sobre a comunidade: a mudança entre os grupos está mais em **quem está presente** ou em **quanto de cada um**?
6. Por que a amostra Mock foi deixada de fora da comparação?

## Se algo der errado

| Problema | O que fazer |
| --- | --- |
| `ls relatorios` não mostra os relatórios | Faça antes a atividade com o Kraken2 (passos 1 a 4) ou use os dados de exemplo do MicrobiomeAnalyst. |
| `metadados.txt` tem menos de 19 amostras | Confira `ls relatorios`: falta o relatório de alguma amostra. Rode `bash classifica.sh` de novo. |
| O site diz que os nomes das amostras não correspondem | Os arquivos não podem ser editados em planilha antes do envio (a planilha pode trocar tabulações por vírgulas). Rode `bash monta_tabelas.sh` de novo e envie os arquivos originais. |
| O site diz que o formato da taxonomia não foi reconhecido | Confira se a opção de **taxonomia em arquivo separado** foi marcada e se o arquivo enviado nesse campo é o `tabela_taxonomia.txt`. |
| O site está lento ou fora do ar | Tente mais tarde. O MicrobiomeAnalyst é um serviço público e pode ficar sobrecarregado nos horários de pico. |
