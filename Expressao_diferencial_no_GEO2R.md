# Expressão diferencial com dados públicos no GEO2R: roteiro passo a passo

## O que vamos fazer

Os glicocorticoides (como a **dexametasona**) são os anti-inflamatórios mais usados no tratamento da asma. Vamos descobrir **que genes a dexametasona liga e desliga nas células do músculo liso das vias aéreas humanas**, usando dados públicos do banco **GEO** (aula: "Dados públicos").

O estudo (Himes et al., 2014) tratou células de **4 doadores** com dexametasona, com albuterol (um broncodilatador), com os dois ou com nada (controle), e sequenciou o RNA das 16 amostras. Vamos comparar apenas **controle × dexametasona** (4 amostras de cada).

Usaremos o **GEO2R**, uma ferramenta do próprio GEO que roda o **DESeq2** no servidor do NCBI a partir das contagens já calculadas. Não é preciso instalar nada: tudo acontece no navegador. A última parte, opcional, usa o terminal para conferir a tabela completa.

O foco não é rodar a análise, que leva um clique, e sim **ler e interpretar** cada gráfico que ela produz, os mesmos vistos em aula.

## 1. Conhecer o estudo no GEO

1. Abra [ncbi.nlm.nih.gov/geo](https://www.ncbi.nlm.nih.gov/geo/) e busque o número de acesso **GSE52778**.
2. Leia o título, o resumo (*Summary*) e o desenho do experimento (*Overall design*).
3. No fim da página, veja a lista de amostras (*Samples*). Cada uma tem um número de acesso `GSM...` (lembre da aula 1: números de acesso). Clique em uma delas para ver as informações da amostra: doador (*cell line*) e tratamento.
4. Na página do GSE, procure o link para o **SRA** (*Relations* ou *SRA Run Selector*): é lá que ficam as reads brutas (`SRR...`).

## 2. Montar a comparação no GEO2R

1. Volte à página do **GSE52778** e clique no botão **Analyze with GEO2R** (logo acima da lista de arquivos).
2. Na página do GEO2R, clique em **Define groups** e crie dois grupos, digitando o nome e apertando **Enter**:
   - `controle`
   - `dexametasona`
3. Na lista de amostras, use as colunas de descrição (título, tratamento) para identificar o tratamento de cada uma. Selecione as **4 amostras sem tratamento** (*untreated*) clicando nas linhas, e depois clique no nome do grupo `controle` para atribuí-las.
4. Selecione as **4 amostras tratadas só com dexametasona** e atribua ao grupo `dexametasona`.
5. **Não coloque em nenhum grupo** as amostras com albuterol (sozinho ou junto com dexametasona). Amostras sem grupo ficam fora da análise.
6. Confira: 4 amostras em cada grupo e 8 amostras sem grupo.
7. Mantenha as opções padrão (aba **Options**) e clique em **Analyze**. A análise leva menos de um minuto.

## 3. Ler os gráficos

Os resultados aparecem em abas com diferentes gráficos e com uma tabela dos genes mais significativos. Abra cada gráfico, observe e anote.

| Gráfico no GEO2R | O que é (aula) | O que observar |
| --- | --- | --- |
| **Histograma de p ajustados** | "Testar milhares de genes ao mesmo tempo" | Há um pico perto de zero? (O GEO2R mostra o p **ajustado**: os genes sem efeito se acumulam perto de 1, em vez de formar o chão plano do histograma de p visto em aula.) |
| **Volcano plot** | "Volcano plot" | Quantos genes ficam nos cantos de cima? Há mais genes induzidos ou reprimidos? |
| **Mean difference (MD)** | "MA plot" | Como os genes pouco expressos se espalham no eixo y? |
| **UMAP** | "Primeiro gráfico: as amostras se agrupam?" | As amostras se agrupam pelo **tratamento** ou pelo **doador**? |
| **Boxplot** | "Por que normalizar?" | As distribuições das amostras ficam parecidas depois da normalização? |

Para descobrir de qual doador é cada ponto do UMAP, passe o mouse sobre ele e compare o código `GSM` com a lista de amostras.

## 4. Ler a tabela de resultados

A tabela embaixo dos gráficos mostra os 250 genes mais significativos, com colunas equivalentes às do DESeq2 vistas em aula (`padj`, `pvalue`, `log2FoldChange`, `baseMean`, nome do gene etc.).

1. Ordene a tabela clicando no nome da coluna `log2FoldChange`.
2. **Confira o sentido da comparação.** O gene **FKBP5** é um alvo clássico dos glicocorticoides e sempre aumenta com a dexametasona. Procure o FKBP5 na tabela (Ctrl + F no navegador). Se o `log2FoldChange` dele for positivo, então valores positivos significam "maior na dexametasona". Se for negativo, a comparação está invertida e o sinal de todos os genes deve ser lido ao contrário.
3. Clique na linha de um gene: aparece um gráfico com a expressão (em TPM) daquele gene em cada amostra.
4. Procure os genes abaixo, todos descritos como respondendo a glicocorticoides, e anote o `log2FoldChange` e o `padj` de cada um. Veja também o gráfico por amostra:

| Gene | Função |
| --- | --- |
| `FKBP5` | regula o próprio receptor de glicocorticoide (retroalimentação negativa) |
| `TSC22D3` | proteína GILZ, anti-inflamatória |
| `PER1` | gene do relógio circadiano |
| `KLF15` | fator de transcrição; metabolismo |
| `ZBTB16` | fator de transcrição |
| `CRISPLD2` | gene destacado no artigo; modula citocinas |

5. Escolha um gene da tabela que você não conhece e procure a função dele no [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene) ou no [UniProt](https://www.uniprot.org).

## 5. O código por trás da análise

Abra a aba **R script**. É o código em R que o GEO2R executou. Procure as linhas com `DESeqDataSetFromMatrix`, `DESeq` e `results`, e compare com o código mostrado na aula (lâmina "Na prática: DESeq2"). Procure também a linha com `design =`.

## 6. Opcional: a tabela completa no terminal

1. No GEO2R, clique em **Download full table** e salve o arquivo na sua Pasta pessoal, dentro de uma pasta chamada `geo2r`.
2. No terminal, entre na pasta e veja o número de cada coluna (troque `NOME_DO_ARQUIVO` pelo nome do arquivo baixado, que aparece com o `ls`):

```bash
cd ~/geo2r
ls
head -n 1 NOME_DO_ARQUIVO | tr '\t' '\n' | cat -n
```

3. Anote o número da coluna `padj` e o da coluna `log2FoldChange`. Conte quantos genes têm p ajustado < 0,05, trocando `N` pelo número da coluna `padj`:

```bash
cd ~/geo2r
awk -F'\t' 'NR > 1 && $N != "" && $N < 0.05' NOME_DO_ARQUIVO | wc -l
```

O `$N != ""` deixa de fora os genes sem `padj`, que o DESeq2 filtra por terem contagens muito baixas (aula: "Lendo a tabela de resultados").

## Perguntas

1. O histograma de p ajustados tem um pico perto de zero? O que isso indica sobre o efeito da dexametasona?
2. No UMAP, as amostras se agruparam pelo tratamento ou pelo doador? O que isso diz sobre a variação entre pessoas diferentes em comparação com o efeito do remédio?
3. Cada doador forneceu uma amostra controle e uma tratada (desenho **pareado**). Olhe a linha `design =` do R script: o doador foi incluído no modelo? Como você escreveria o `design` para incluí-lo? Por que isso poderia aumentar o número de genes detectados? (Lembre da lâmina sobre efeito de lote.)
4. Os 6 genes da tabela do passo 4.4 aumentaram, diminuíram ou não mudaram? Algum deles é induzido em todos os 4 doadores, no gráfico por amostra?
5. Há mais genes induzidos ou reprimidos pela dexametasona no volcano plot? Os glicocorticoides agem ligando-se a um receptor que é um fator de transcrição: isso ajuda a explicar o padrão?
6. Um gene muito induzido nestas células em cultura é necessariamente importante para o efeito do remédio no pulmão de um paciente? Cite dois motivos para cautela (lâmina "Antes de concluir").
7. O GEO2R usou contagens calculadas pelo NCBI a partir das reads do SRA. Por que é importante que os autores tenham depositado as reads brutas, e não só a lista final de genes?

## Se algo der errado

| Problema | O que fazer |
| --- | --- |
| Não aparece o botão **Analyze with GEO2R** | Confira se o acesso digitado é exatamente **GSE52778**. |
| O GEO2R pede para escolher uma plataforma | Escolha a plataforma de sequenciamento listada (só há uma para este estudo). |
| O botão **Analyze** não fica ativo | É preciso ter pelo menos uma amostra em cada um dos dois grupos. Confira a atribuição dos grupos. |
| A análise dá erro ou não termina | Recarregue a página e monte os grupos de novo. Se persistir, tente mais tarde: o GEO2R roda nos servidores do NCBI. |
| O FKBP5 aparece com `log2FoldChange` negativo | A comparação está na ordem inversa. Leia todos os sinais ao contrário (passo 4.2). |
| No terminal, o `awk` devolve `0` | Confira o número da coluna `padj` com o comando `head -n 1 ... | tr ...` e se trocou o `N` pelo número. |

## Referência

Himes BE et al., 2014. RNA-Seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid responsive gene that modulates cytokine function in airway smooth muscle cells. PLoS One 9: e99625.
