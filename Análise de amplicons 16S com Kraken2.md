# Perfil taxonômico de amplicons 16S com Kraken2: roteiro passo a passo

## O que vamos fazer

Vamos descobrir **quais bactérias vivem no intestino de um camundongo** e como essa comunidade muda com a idade, usando reads de sequenciamento do gene **16S rRNA** (região V4, Illumina MiSeq 2×250 pb) amplificado por PCR.

Os dados são um conjunto de exemplo. Contém amostras de fezes de uma mesma fêmea de camundongo (F3) coletadas em dois períodos:

| Amostras | Período |
| --- | --- |
| `F3D0` a `F3D9` | primeiros 10 dias após o desmame |
| `F3D141` a `F3D150` | cerca de 5 meses depois |
| `Mock` | **comunidade artificial**: DNA misturado de cerca de 20 bactérias conhecidas |

A amostra "**Mock"** é o nosso controle: como sabemos o que foi colocado nela, podemos conferir se a análise está correta.

Para classificar as reads vamos usar o **Kraken2** com um banco de referência de 16S (SILVA). O Kraken2 quebra cada read em pequenos pedaços (*k-mers*) e procura esses pedaços no banco, atribuindo a read ao táxon mais provável.

Todos os arquivos ficam na pasta `metagenomica`, na sua pasta pessoal.

## 1. Instalar os programas

Se o micromamba ainda não estiver instalado, siga primeiro a seção 1 do roteiro de filogenômica (instalação do micromamba e `source ~/.bashrc`).

Abra o terminal (pode ser usado o atalho **Ctrl + Alt + T**) e crie um ambiente micromamba novo, chamado `metagenomica`, com os programas desta atividade:

> **Comandos em várias linhas:** comandos longos foram quebrados em linhas curtas. Uma barra invertida `\` no fim da linha indica que o comando continua na linha de baixo. Copie e cole o bloco inteiro de uma vez.

```bash
micromamba create -n metagenomica \
  -c conda-forge -c bioconda -y \
  kraken2 seqkit
```

Ative o ambiente e crie a pasta da atividade:

```bash
micromamba activate metagenomica
mkdir ~/metagenomica
cd ~/metagenomica
```

> **Importante:** sempre que abrir um terminal novo, rode `micromamba activate metagenomica` e `cd ~/metagenomica` antes de continuar.

## 2. Baixar os dados e o banco de referência

### 2.1 Reads de sequenciamento

O comando `wget` baixa um arquivo da internet direto para a pasta em que você está. Depois, o `unzip` extrai o `.zip` e o `m``v` procura todos os arquivos `.fastq` em qualquer pasta (com o uso do `*`) no local atual e os move para a pasta `reads`:

```bash
cd ~/metagenomica
wget https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip
unzip -o miseqsopdata.zip
mkdir reads
mv */*.fastq reads/
ls reads
```

Devem aparecer 40 arquivos: um par (`R1` e `R2`) para cada uma das 20 amostras, por exemplo `F3D0_S188_L001_R1_001.fastq` e `F3D0_S188_L001_R2_001.fastq`.

### 2.2 Banco de referência do 16S (SILVA)

O banco já vem pronto para o Kraken2. O download pode levar alguns minutos:

```bash
cd ~/metagenomica
wget https://genome-idx.s3.amazonaws.com/kraken/16S_Silva138_20200326.tgz
mkdir banco
tar -xzf 16S_Silva138_20200326.tgz -C banco
```

Um banco do Kraken2 é uma pasta com arquivos que terminam em `.k2d`. Como não sabemos exatamente em que subpasta eles ficaram depois da extração, pedimos para o `find` encontrar o arquivo `hash.k2d` e guardamos o nome da pasta na variável `DB`:

```bash
cd ~/metagenomica
DB=$(dirname $(find banco -name hash.k2d))
echo $DB
```

Deve aparecer o caminho da pasta do banco (algo como `banco/16S_SILVA138_k2db`). Essa variável com o caminho da DB será usada mais tarde.

## 3. Conhecendo os dados

Antes de analisar, olhe os arquivos. Lembre-se do formato FASTQ (aula 2): 4 linhas por read, com a qualidade de cada base na quarta linha.

```bash
cd ~/metagenomica
head -n 4 reads/F3D0_S188_L001_R1_001.fastq
seqkit stats reads/*.fastq
```

O `seqkit stats` mostra, para cada arquivo, o número de reads (`num_seqs`) e o tamanho médio e máximo delas (`avg_len`, `max_len`).

## 4. Classificação das reads

### 4.1 Entendendo com um exemplo: uma amostra

Antes de rodar tudo, faça o processo uma vez à mão com a amostra `F3D0`. Se você abriu um terminal novo, defina a variável `DB` de novo (último comando do passo 2.2).

```bash
cd ~/metagenomica
kraken2 --db $DB --paired --threads 2 \
  --use-mpa-style --report teste.txt \
  --output /dev/null \
  reads/F3D0_S188_L001_R1_001.fastq \
  reads/F3D0_S188_L001_R2_001.fastq
```

O que cada opção faz:

- `--db $DB`: a pasta do banco de referência (uso de variável para indicar o local do banco não é obrigatório).
- `--paired`: os dois arquivos são as reads R1 e R2 de um mesmo fragmento (*paired-end*).
- `--threads 2`: usa 2 processadores.
- `--report teste.txt` e `--use-mpa-style`: escreve um relatório com quantas reads caíram em cada táxon, com a linhagem completa (domínio, filo, classe, ordem, família, gênero).
- `--output /dev/null`: descarta o resultado read por read, que não vamos usar.

Ao terminar, o Kraken2 escreve na tela quantas sequências foram classificadas e qual a porcentagem.

Agora veja o relatório resultante:

```bash
cd ~/metagenomica
cat teste.txt
```

Cada linha tem uma linhagem e o número de reads, por exemplo:

```
d__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Muribaculaceae|g__Muribaculaceae    2531
```

Os prefixos indicam o nível taxonômico: `d__` domínio, `p__` filo, `c__` classe, `o__` ordem, `f__` família, `g__` gênero. Esse é um padrão usado em muitos programas de metagenômica.

Opcional: Checar os 10 gêneros mais abundantes usando o `awk`. Ele lê o arquivo linha por linha e, para as linhas que **terminam** em um gênero (`g__`), imprime o número de reads e só o último nome da linhagem. O `sort -nr` ordena do maior para o menor número:

```bash
cd ~/metagenomica
awk -F'\t' '$1 ~ /g__[^|]*$/ {n = split($1, nome, "|"); print $2 "\t" nome[n]}' \
  teste.txt | sort -nr | head -n 10
```

### 4.2 Todas as amostras

Agora repetimos a classificação para as 20 amostras com um loop `for`, dentro de um **script**. Dois detalhes:

- `${r1/_R1_/_R2_}` troca `_R1_` por `_R2_` no nome do arquivo que está dentro da variável`r1`, para achar o par de cada R1;
- `${amostra%%_*}` apaga tudo a partir do primeiro `_`, deixando só o nome da amostra (`F3D0_S188_L001_R1_001.fastq` vira `F3D0`).

1. Abra o **Editor de Texto** do Ubuntu.
2. Copie e cole todo o conteúdo abaixo.

```bash
#!/bin/bash
# Classifica as reads de cada amostra com o Kraken2.
# Relatorios vao para a pasta relatorios/ e
# o resumo de cada corrida vai para a pasta logs/.

DB=$(dirname $(find banco -name hash.k2d))
mkdir -p relatorios logs

for r1 in reads/*_R1_001.fastq; do
  r2=${r1/_R1_/_R2_}
  amostra=$(basename $r1)
  amostra=${amostra%%_*}
  echo "Classificando $amostra..."
  kraken2 --db $DB --paired --threads 2 \
    --use-mpa-style --report relatorios/$amostra.txt \
    --output /dev/null \
    $r1 $r2 2> logs/$amostra.log
done

echo "Pronto!"
```

3. Salve o arquivo na pasta **`metagenomica`** com o nome **`classifica.sh`**.
4. No terminal, execute o script:

```bash
cd ~/metagenomica
bash classifica.sh
```

O `2> logs/$amostra.log` guarda o resumo que o Kraken2 escreveria na tela. Para ver a porcentagem de reads classificadas de todas as amostras:

```bash
cd ~/metagenomica
ls relatorios
grep "classified" logs/*.log
```

## 5. O controle: a comunidade Mock

A amostra Mock foi feita misturando DNA de bactérias conhecidas (comunidade de referência do *Human Microbiome Project*). Entre elas estão espécies destes gêneros:

| Gêneros esperados na Mock |
| --- |
| *Acinetobacter*, *Actinomyces*, *Bacillus*, *Bacteroides*, *Clostridium*, *Deinococcus* |
| *Enterococcus*, *Escherichia*, *Helicobacter*, *Lactobacillus*, *Listeria*, *Neisseria* |
| *Propionibacterium*, *Pseudomonas*, *Rhodobacter*, *Staphylococcus*, *Streptococcus* |

Veja os gêneros encontrados na Mock (o mesmo comando do passo 4.1, agora com `head -n 30`):

```bash
cd ~/metagenomica
awk -F'\t' '$1 ~ /g__[^|]*$/ {n = split($1, nome, "|"); print $2 "\t" nome[n]}' \
  relatorios/Mock.txt | sort -nr | head -n 30
```

Compare com a tabela. Atenção: alguns nomes mudaram desde a montagem da comunidade e o banco SILVA pode usar o nome novo ou um nome composto (por exemplo, `Escherichia-Shigella`). Quando um gênero esperado "faltar", procure o nome da espécie no [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) para ver se ele foi reclassificado.

## 6. Diversidade alfa

### 6.1 Subamostragem (rarefação)

As amostras têm números diferentes de reads, e quanto mais reads, mais gêneros pode ser possível observar por consequência de um maior "esforço amostral" e não por um motivo biológico. Para comparar a diversidade de forma justa, vamos **sortear o mesmo número de reads de cada amostra** e classificar de novo. É a rarefação vista em aula.

Primeiro, descubra qual amostra tem menos reads. O `-T` deixa a saída do `seqkit stats` em colunas, e o `sort` ordena pela quarta coluna (`num_seqs`), do menor para o maior:

```bash
cd ~/metagenomica
seqkit stats -T reads/*_R1_001.fastq \
  | sort -t$'\t' -k4,4n
```

Anote o menor `num_seqs` entre as amostras `F3D` e escolha um número redondo um pouco abaixo dele: esse será o `N`, o número de pares de reads sorteados por amostra.

Abra o **Editor de Texto**, cole o conteúdo abaixo, troque o `2000` da linha `N=2000` pelo número que você escolheu e salve na pasta `metagenomica` com o nome **`subamostra.sh`**. O `seqkit sample` sorteia as reads; usamos a mesma **seed** (`-s 11`) no R1 e no R2 para que sejam sorteados os mesmos pares nos dois arquivos.

```bash
#!/bin/bash
# Sorteia N pares de reads de cada amostra (rarefacao)
# e classifica de novo com o Kraken2.

N=2000   # pares de reads por amostra (ajuste no passo 6.1)

DB=$(dirname $(find banco -name hash.k2d))
mkdir -p reads_sub relatorios_sub logs_sub

for r1 in reads/*_R1_001.fastq; do
  r2=${r1/_R1_/_R2_}
  amostra=$(basename $r1)
  amostra=${amostra%%_*}
  echo "Subamostrando e classificando $amostra..."
  seqkit sample -n $N -s 11 $r1 \
    -o reads_sub/${amostra}_R1.fastq 2> /dev/null
  seqkit sample -n $N -s 11 $r2 \
    -o reads_sub/${amostra}_R2.fastq 2> /dev/null
  kraken2 --db $DB --paired --threads 2 \
    --use-mpa-style --report relatorios_sub/$amostra.txt \
    --output /dev/null \
    reads_sub/${amostra}_R1.fastq reads_sub/${amostra}_R2.fastq \
    2> logs_sub/$amostra.log
done

echo "Pronto!"
```

Execute e confira se todas as amostras ficaram com  `N` reads na coluna `n``u``m``_seqs`, e se o R1 e o R2 de cada amostra têm o mesmo número:

```bash
cd ~/metagenomica
bash subamostra.sh
seqkit stats reads_sub/*.fastq
```

### 6.2 Riqueza e índice de Shannon

Riqueza e índice de Shannon são duas medidas de diversidade alfa que respondem a perguntas diferentes sobre a mesma amostra:

- **Riqueza:** quantos gêneros diferentes existem na amostra. Cada gênero conta uma vez só, tenha ele 1 read ou 10 mil.
- **Índice de Shannon:** leva em conta quantos gêneros existem **e** como as reads se distribuem entre eles (a equitabilidade). Ele aumenta quando há mais gêneros e quando nenhum gênero domina a amostra.

Consequência prática: a riqueza depende muito dos gêneros raros, que só aparecem quando há muitas reads. Por isso ela é mais afetada pela profundidade de sequenciamento do que o Shannon (lembre-se do passo 6.1).

Para o cálculo de diversidade vamos usar o **R** com o pacote **vegan**, um dos pacotes mais usados em ecologia de comunidades. A riqueza vem da função `specnumber()` e o Shannon da função `diversity()`.

Instale o R e o vegan no ambiente da atividade (leva alguns minutos):

```bash
micromamba activate metagenomica
micromamba install -n metagenomica \
  -c conda-forge -y r-base r-vegan
```

Abra o **Editor de Texto**, cole o conteúdo abaixo e salve na pasta `metagenomica` com o nome **`diversidade.R`**. As linhas que começam com `#` são comentários.

```r
# Calcula a riqueza e o indice de Shannon (nivel de genero)
# de cada amostra, usando o pacote vegan.
# Uso: Rscript diversidade.R PASTA   (relatorios ou relatorios_sub)

library(vegan)

pasta <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(pasta)) stop("Indique a pasta. Exemplo: Rscript diversidade.R relatorios")
arquivos <- list.files(pasta, pattern = "\\.txt$", full.names = TRUE)

# 1) Ler os relatorios e montar a tabela amostra x genero
contagens <- list()
for (arquivo in arquivos) {
  amostra <- sub("\\.txt$", "", basename(arquivo))
  tab <- read.delim(arquivo, header = FALSE, quote = "",
                    col.names = c("linhagem", "reads"))
  tab <- tab[grepl("g__[^|]*$", tab$linhagem), ]
  contagens[[amostra]] <- setNames(tab$reads, tab$linhagem)
}
generos <- unique(unlist(lapply(contagens, names)))
comunidade <- matrix(0, nrow = length(contagens), ncol = length(generos),
                     dimnames = list(names(contagens), generos))
for (amostra in names(contagens)) {
  comunidade[amostra, names(contagens[[amostra]])] <- contagens[[amostra]]
}

# 2) Riqueza e Shannon com o vegan
dia <- suppressWarnings(as.numeric(sub("^F3D", "", rownames(comunidade))))
grupo <- ifelse(is.na(dia), "Mock", ifelse(dia <= 9, "Dias_0a9", "Dias_141a150"))
resultado <- data.frame(
  amostra = rownames(comunidade),
  grupo   = grupo,
  riqueza = specnumber(comunidade),
  shannon = round(diversity(comunidade, index = "shannon"), 3)
)
write.table(resultado, paste0("diversidade_", pasta, ".txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(resultado, row.names = FALSE)

# 3) Media de cada grupo e boxplots (sem a Mock)
amostras_f3 <- resultado[resultado$grupo != "Mock", ]
print(aggregate(cbind(riqueza, shannon) ~ grupo, data = amostras_f3, FUN = mean))
pdf(paste0("diversidade_", pasta, ".pdf"), width = 8, height = 4)
par(mfrow = c(1, 2))
boxplot(riqueza ~ grupo, data = amostras_f3, main = "Riqueza", ylab = "Generos")
boxplot(shannon ~ grupo, data = amostras_f3, main = "Shannon", ylab = "Shannon")
dev.off()
```

O que cada bloco faz:

- **Bloco 1:** lê os relatórios do Kraken2, guarda só as linhas que terminam em um gênero (`g__`) e monta a **tabela amostra × gênero** vista na aula, com `0` quando o gênero não apareceu na amostra.
- **Bloco 2:** usa o nome da amostra para definir o grupo (dia até 9, dia 141 em diante ou Mock) e calcula a riqueza e o Shannon com o vegan.
- **Bloco 3:** mostra a média de cada grupo e salva os boxplots em PDF.

Execute o script duas vezes, uma para cada pasta de relatórios:

```bash
cd ~/metagenomica
Rscript diversidade.R relatorios
Rscript diversidade.R relatorios_sub
```

Cada execução cria dois arquivos: uma tabela (`diversidade_relatorios.txt` com todas as reads e `diversidade_relatorios_sub.txt` com as reads subamostradas) e um PDF com os boxplots dos dois grupos (`diversidade_relatorios.pdf` e `diversidade_relatorios_sub.pdf`). Abra os PDFs e compare os resultados com e sem subamostragem.

## 7. Perguntas

1. Quais são os 5 gêneros mais abundantes na amostra `F3D0` e na amostra `F3D150`? A comunidade mudou?
2. Quantos dos gêneros esperados foram encontrados na Mock? Apareceu algum gênero que **não** deveria estar lá? O que isso pode indicar?
3. Compare a riqueza e o Shannon dos dois grupos de amostras. Qual grupo é mais diverso?
4. O número de reads varia entre as amostras (veja o `seqkit stats`). Compare a riqueza e o Shannon calculados com todas as reads e com as reads subamostradas. Qual das duas medidas mudou mais, e por quê?
5. Por que o banco SILVA classifica bem até gênero, mas raramente até espécie, com reads da região V4?

## Se algo der errado

| Problema | O que fazer |
| --- | --- |
| `command not found` | Rode `micromamba activate metagenomica`. |
| `No such file or directory` | Rode `cd ~/metagenomica` e `ls`. Confira se existem as pastas `reads` e `banco`. |
| `ls reads` não mostra os 40 arquivos `.fastq` | Repita o passo 2.1. Se o `wget` falhar, abra o link no navegador, salve o `.zip` na pasta `metagenomica` e continue a partir do `unzip`. |
| `echo $DB` não mostra nada | O banco não foi extraído. Confira com `ls banco` e repita o passo 2.2. |
| Kraken2 diz `database ("...") not found` | A variável `DB` não existe neste terminal. Rode de novo `DB=$(dirname $(find banco -name hash.k2d))`. |
| Algum relatório em `relatorios/` está vazio | Veja o arquivo `.log` da amostra na pasta `logs` (`cat logs/F3D0.log`). |
| `Killed` ou o computador trava durante o Kraken2 | Pouca memória. Feche outros programas e troque `--threads 2` por `--threads 1`. |
