# Seleção natural com dN/dS: roteiro passo a passo

## O que vamos fazer

Quando um gene de duas espécies é comparado, as diferenças no DNA podem ser de dois tipos:

- **não sinônimas (dN):** mudam o aminoácido da proteína;
- **sinônimas (dS):** não mudam o aminoácido (por causa da redundância do código genético).

A razão **dN/dS** (também chamada de ω) indica que tipo de seleção natural agiu sobre o gene:

| dN/dS | Interpretação |
| --- | --- |
| bem menor que 1 | **seleção purificadora**: mudanças na proteína são eliminadas porque prejudicam a função |
| perto de 1 | evolução neutra: mudar ou não o aminoácido dá no mesmo |
| maior que 1 | **seleção positiva**: mudanças na proteína foram favorecidas |

Vamos calcular o dN/dS de vários genes comparando duas leveduras muito próximas: ***Saccharomyces cerevisiae*** e ***Saccharomyces paradoxus***. Usamos espécies próximas porque, entre espécies muito distantes, as mutações sinônimas se acumulam tanto que deixam de ser contadas corretamente (o dS "satura").

**Genes:** as mesmas 8 proteínas da atividade de filogenômica (RPO21, RPB2, RPA190, PRP8, HSH155, MCM2, MCM7, TCP1). Caso não tenha seguido o passo a passo da filogenômica, leia a seção do guia sobre como baixar proteínas no Uniprot.

**Hipótese para testar.** Genes **essenciais** (sem os quais a levedura morre) devem estar sob seleção purificadora mais forte, e portanto ter dN/dS menor, que genes **não essenciais** (Hirsh & Fraser, 2001). As 8 proteínas acima são essenciais (grupo A). Para contrastar, baixe também as 6 proteínas abaixo, que não são essenciais em laboratório (grupo B), do mesmo jeito que as do grupo A:

| Gene (grupo B) | Função |
| --- | --- |
| GAL4 | ativa os genes de uso da galactose |
| GAL80 | reprime o GAL4 quando não há galactose |
| GAL7 | enzima do uso da galactose |
| GAL10 | enzima do uso da galactose |
| STE2 | receptor do feroônio de acasalamento α |
| STE3 | receptor do feroônio de acasalamento a |

O script calcula todos de uma vez. Na hora de montar a tabela de resultados, indique a qual grupo (A ou B) cada gene pertence.

Todos os arquivos ficam na pasta `dnds`, na sua Pasta pessoal.

## 1. Instalar os programas

Se o micromamba ainda não estiver instalado, siga primeiro a seção 1 do roteiro de filogenômica (instalação do micromamba e `source ~/.bashrc`).

Abra o terminal (**Ctrl + Alt + T**) e crie um ambiente novo, chamado `dnds`, com os programas desta atividade. Isso leva alguns minutos:

> **Comandos em várias linhas:** comandos longos foram quebrados em linhas curtas. Uma barra invertida `\` no fim da linha indica que o comando continua na linha de baixo. Copie e cole o bloco inteiro de uma vez.

```bash
micromamba create -n dnds \
  -c conda-forge -c bioconda -y \
  blast seqkit mafft pal2nal paml
```

O que cada um faz: **BLAST** encontra o gene em cada espécie; **SeqKit** separa sequências; **MAFFT** alinha as proteínas; **PAL2NAL** transforma o alinhamento de proteínas em alinhamento de códons (DNA); **PAML** (programa `codeml`) calcula dN e dS.

Ative o ambiente e crie a pasta da atividade:

```bash
micromamba activate dnds
mkdir ~/dnds
cd ~/dnds
```

> **Importante:** sempre que abrir um terminal novo, rode `micromamba activate dnds` e `cd ~/dnds` antes de continuar.

## 2. Baixar os arquivos

Para calcular dN/dS precisamos, para cada espécie, das **proteínas** (para encontrar e alinhar os genes) e das **sequências de DNA codificante** (as CDS, onde contamos as mutações sinônimas e não sinônimas).

### 2.1 Genomas das duas espécies (NCBI)

Repita para *Saccharomyces cerevisiae* e para *Saccharomyces paradoxus*:

1. Abra [ncbi.nlm.nih.gov/datasets/genome](https://www.ncbi.nlm.nih.gov/datasets/genome/) e digite o nome da espécie.
2. Clique no genoma marcado como **reference genome** (código começando com **GCF\_**).
3. Clique em **Download** e marque **duas opções**: **Protein (FASTA)** e **Genomic coding sequences (FASTA)**. Baixe o `.zip`.
4. Extraia o `.zip` (botão direito > Extrair aqui) e entre em `ncbi_dataset` > `data` > `GCF_...`. Lá estão dois arquivos: `protein.faa` (proteínas) e `cds_from_genomic.fna` (CDS).
5. Copie os dois para a pasta `dnds` e renomeie com **o mesmo nome** (espécie + código GCF, sem espaços), mudando só a terminação:

| Arquivo original | Novo nome (exemplo) |
| --- | --- |
| `protein.faa` | `Saccharomyces_cerevisiae_GCF_000146045.2.faa` |
| `cds_from_genomic.fna` | `Saccharomyces_cerevisiae_GCF_000146045.2.fna` |

Use o código GCF que aparece no site para cada espécie.

### 2.2 Proteínas de referência (UniProt)

Use as mesmas proteínas da filogenômica. Se já tem os arquivos `RPO21.fasta`, `RPB2.fasta` etc., copie-os para a pasta `dnds`. Se não tem, baixe cada uma no [UniProt](https://www.uniprot.org): busque, por exemplo, `RPO21 Saccharomyces cerevisiae`, escolha o resultado revisado de *S. cerevisiae* (S288c), clique em **Download** > **FASTA** e salve como `RPO21.fasta`.

Para analisar mais genes, basta colocar mais arquivos `.fasta` na pasta.

### 2.3 Conferir

```bash
cd ~/dnds
ls
```

Devem aparecer: 2 arquivos `.faa`, 2 arquivos `.fna` (com os mesmos nomes dos `.faa`) e os arquivos `.fasta` das proteínas de referência.

## 3. Entender os passos e criar os arquivos

### 3.1 O que acontece com cada proteína

Para cada proteína de referência (cada arquivo `.fasta`), o script vai fazer:

- **Passo A — BLAST:** procura a proteína em cada espécie e guarda o código do melhor resultado (por exemplo, `NP_010141.1`).
- **Passo B — proteína:** copia a sequência dessa proteína para o arquivo `par.faa`, com o nome da espécie.
- **Passo C — CDS:** procura no arquivo `.fna` a sequência de DNA que codifica essa mesma proteína (o cabeçalho de cada CDS tem o texto `protein_id=` seguido do código da proteína), tira o códon de parada (os 3 últimos nucleotídeos) e copia para `par.fna`.
- **Passo D — alinhamento:** o MAFFT alinha as duas proteínas; o PAL2NAL usa esse alinhamento para alinhar o DNA códon por códon, e salva no formato do PAML (`par.phy`).
- **Passo E — dN/dS:** o `codeml` compara as duas sequências e calcula dN, dS e dN/dS. O resultado de cada proteína vai para o arquivo `resultados.txt`.

Os arquivos `par.*` são temporários: são apagados e refeitos a cada proteína.

### 3.2 Criar o arquivo de configuração do codeml

O `codeml` lê as opções de um arquivo de texto. Abra o **Editor de Texto**, cole o conteúdo abaixo e salve na pasta `dnds` com o nome **`codeml.ctl`**:

```
seqfile = par.phy
outfile = par.mlc
runmode = -2
seqtype = 1
CodonFreq = 2
model = 0
NSsites = 0
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.4
cleandata = 1
```

As linhas mais importantes: `seqfile` e `outfile` são os arquivos de entrada e saída; `runmode = -2` compara as sequências aos pares; `seqtype = 1` diz que são códons; `icode = 0` usa o código genético padrão.

### 3.3 Criar o script

Abra outro arquivo no **Editor de Texto**, cole o conteúdo abaixo e salve na pasta `dnds` com o nome **`calcula_dnds.sh`**. A primeira linha (`#!/bin/bash`, o *shebang*) avisa que o arquivo deve ser lido pelo bash. As linhas que começam com `#` são comentários e mostram onde está cada passo.

O script não tem nenhum nome de espécie ou de gene escrito: ele usa todos os arquivos `.faa` (espécies) e `.fasta` (proteínas) que estiverem na pasta.

```bash
#!/bin/bash
# Calcula dN/dS entre as duas especies (.faa e .fna)
# para cada proteina de referencia (.fasta).

echo "Resultados do dN/dS" > resultados.txt

# Bancos do BLAST, um por especie
for genoma in *.faa; do
  makeblastdb -in $genoma -dbtype prot > /dev/null
done

for proteina in *.fasta; do
  rm -f par.faa par.fna par.aln par.phy par.mlc

  for genoma in *.faa; do
    nome=${genoma%.faa}
    # Passo A: codigo do melhor hit nesta especie
    blastp -query $proteina -db $genoma \
      -evalue 1e-10 -max_target_seqs 1 -max_hsps 1 \
      -outfmt "6 sseqid" > melhor_hit.txt
    id=$(cat melhor_hit.txt)
    # Passo B: sequencia da proteina encontrada
    seqkit grep -f melhor_hit.txt $genoma \
      | sed "s/^>.*/>$nome/" >> par.faa
    # Passo C: CDS da mesma proteina, sem o codon de parada
    seqkit grep -n -r -p "protein_id=$id\]" $nome.fna \
      | seqkit subseq -r 1:-4 \
      | sed "s/^>.*/>$nome/" >> par.fna
  done

  # Passo D: alinhar proteinas e converter em codons
  mafft --auto --quiet par.faa > par.aln
  pal2nal.pl par.aln par.fna -output paml -nogap \
    > par.phy

  # Passo E: calcular dN/dS e guardar o resultado
  codeml codeml.ctl > /dev/null
  echo "${proteina%.fasta} $(grep 'dN/dS=' par.mlc)" \
    >> resultados.txt
done

echo "Pronto!"
```

## 4. Rodar e ler os resultados

Confira se os dois arquivos que você criou estão na pasta e execute o script:

```bash
cd ~/dnds
ls codeml.ctl calcula_dnds.sh
bash calcula_dnds.sh
```

Quando aparecer `Pronto!`, veja os resultados:

```bash
cd ~/dnds
cat resultados.txt
```

Cada linha traz o nome da proteína e o resultado do `codeml`, mais ou menos assim:

```
RPO21 t= 0.1234  S= 1500.0  N= 4500.0  dN/dS= 0.0150  dN = 0.0030  dS = 0.2000
```

| Campo | Significado |
| --- | --- |
| `S` e `N` | número de sítios sinônimos e não sinônimos no gene |
| `dN` | taxa de substituições não sinônimas |
| `dS` | taxa de substituições sinônimas |
| `dN/dS` | a razão ω (o valor principal desta atividade) |

(Os números acima são só um exemplo do formato.)

Monte uma tabela com o nome de cada proteína, o grupo (A ou B), dN, dS e dN/dS, e faça um gráfico com o dN/dS de cada uma, separando os dois grupos (em uma planilha, por exemplo).

**Confira o dS:** se o dS de algum gene for muito alto (acima de 2), as mutações sinônimas provavelmente saturaram e o dN/dS desse gene não é confiável.

## Se algo der errado

| Problema | O que fazer |
| --- | --- |
| `command not found` | Rode `micromamba activate dnds`. |
| `No such file or directory` | Rode `cd ~/dnds` e `ls`. Confira os nomes: cada `.faa` precisa de um `.fna` com exatamente o mesmo nome. |
| Uma linha de `resultados.txt` tem só o nome da proteína, sem números | O BLAST não achou a proteína em uma das espécies, ou a CDS não foi encontrada. Veja `cat par.faa` e `cat par.fna` logo depois de rodar: os dois devem ter 2 sequências. |
| Todas as linhas sem números | Confira se o arquivo `codeml.ctl` está na pasta e se foi salvo exatamente com esse nome. |
| Mensagens do PAL2NAL sobre sequências que não correspondem | Confira se o `.fna` e o `.faa` de cada espécie vieram do mesmo download (mesmo código GCF). |
