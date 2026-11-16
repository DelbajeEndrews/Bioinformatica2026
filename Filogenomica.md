# Filogenômica com BLAST: roteiro passo a passo

## 1. Instalação do micromamba e dos programas

### Instalação do micromamba

Abra o terminal (**Ctrl + Alt + T**) e cole o comando (colar no terminal é **Ctrl + Shift + V**):

```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

O instalador faz quatro perguntas. Responda assim:

| Pergunta | O que fazer |
| --- | --- |
| `Micromamba binary folder? [~/.local/bin]` | pressionar **Enter** |
| `Init shell (bash)? [Y/n]` | digitar **y** e pressionar **Enter** |
| `Configure conda-forge? [Y/n]` | digitar **y** e pressionar **Enter** |
| `Prefix location? [~/micromamba]` | pressionar **Enter** |

Depois, recarregue as configurações do terminal:

```bash
source ~/.bashrc
```

### Instalação dos programas dentro de um *environment*

O *environment* (ambiente) é uma caixa onde ficam os programas desta atividade. Vamos criá-lo com o nome `filogenomica` e instalar o BLAST, o SeqKit, o MAFFT e o IQ-TREE:

> **Comandos em várias linhas:** comandos longos foram quebrados em linhas curtas. Uma barra invertida `\` no fim da linha indica que o comando continua na linha de baixo. Copie e cole o bloco inteiro de uma vez.

```bash
micromamba create -n filogenomica \
  -c conda-forge -c bioconda -y \
  blast seqkit mafft iqtree
```

Ative o ambiente e crie a pasta onde todos os arquivos vão ficar:

```bash
micromamba activate filogenomica
mkdir ~/filogenomica
cd ~/filogenomica
```

> **Importante:** sempre que abrir um terminal novo, rode `micromamba activate filogenomica` e `cd ~/filogenomica` antes de continuar. Lembre-se: geralmente o diretório "\~" se refere ao local "/home/usuário" (/home/aluno).

## 2. Espécies e download das proteínas

A lista abaixo tem 8 leveduras. *S. cerevisiae* entra sempre, porque é a espécie de referência. Escolham de **2 a 4** espécies do meio da tabela e incluam o **grupo externo** (*S. pombe*), que serve para enraizar a árvore. Baixe de 6 a 10 proteomas (.faa) das espécies escolhidas.

| Espécie | Papel |
| --- | --- |
| *Saccharomyces cerevisiae* | referência (obrigatória) |
| *Saccharomyces paradoxus* |  |
| *Nakaseomyces glabratus* (antiga *Candida glabrata*) |  |
| *Zygosaccharomyces rouxii* |  |
| *Kluyveromyces lactis* |  |
| *Lachancea thermotolerans* |  |
| *Eremothecium gossypii* |  |
| *Schizosaccharomyces pombe* | grupo externo (obrigatório) |

### Como baixar o FASTA de proteínas (.faa) no NCBI

Repita para cada espécie escolhida, incluindo *S. cerevisiae* e *S. pombe*:

1. Abra [ncbi.nlm.nih.gov/datasets/genome](https://www.ncbi.nlm.nih.gov/datasets/genome/) e digite o nome da espécie.
2. Clique no genoma marcado como **reference genome** (código começando com **GCF\_**).
3. Clique em **Download**, marque **somente Protein (FASTA)** e baixe o arquivo `.zip`.
4. Na pasta Downloads (ou Transferências), clique com o botão direito no `.zip` e escolha **Extrair aqui**.
5. Entre nas pastas `ncbi_dataset` > `data` > `GCF_...` e encontre o arquivo **`protein.faa`**. Renomeie o arquivo FASTA com o nome da espécie + o GCF (ou GCA) accession code.
6. Copie esse arquivo renomeado para a sua pasta **`filogenomica`**.

No fim, confira no terminal:

```bash
cd ~/filogenomica
ls
```

Devem aparecer os arquivos `.faa` de todas as espécies escolhidas.

## 3. Proteínas de referência

Vamos procurar 8 proteínas de *S. cerevisiae* em todos os genomas. Elas foram escolhidas porque: (1) existem em todas as leveduras; (2) aparecem em **uma única cópia** por genoma, o que evita confundir genes parecidos (parálogos); (3) têm funções essenciais e evoluem devagar, sendo fáceis de alinhar; e (4) são longas, então carregam bastante informação para a árvore.

Obs: para montar uma lista semelhante para outro conjunto de espécies devem ser feitas buscas na literatura científica + auxílio de IAs.

| Gene | Função |
| --- | --- |
| RPO21 | maior subunidade da RNA polimerase II |
| RPB2 | segunda subunidade da RNA polimerase II |
| RPA190 | maior subunidade da RNA polimerase I |
| PRP8 | núcleo do spliceossomo |
| HSH155 | componente do spliceossomo |
| MCM2 | helicase de replicação do DNA |
| MCM7 | helicase de replicação do DNA |
| TCP1 | chaperonina (ajuda outras proteínas a se dobrarem) |

### Como baixar cada proteína no UniProt

Repita para os 8 genes:

1. Abra [uniprot.org](https://www.uniprot.org) e busque o nome do gene seguido da espécie, por exemplo: `RPO21 Saccharomyces cerevisiae`.
   1. Obs: experimente fazer a mesma busca no NCBI como alternativa.
2. Escolha o resultado **revisado** (ícone com folha dourada) do organismo *Saccharomyces cerevisiae (strain ATCC 204508 / S288c)*.
3. Clique em **Download**, escolha o formato **FASTA** e baixe.
4. Mova o arquivo para a pasta `filogenomica` e renomeie com o nome do gene e a terminação `.fasta`: `RPO21.fasta`, `RPB2.fasta`, `RPA190.fasta`, `PRP8.fasta`, `HSH155.fasta`, `MCM2.fasta`, `MCM7.fasta`, `TCP1.fasta`.

Repare na diferença: os genomas terminam em **`.faa`** (passo 2) e as proteínas de referência em **`.fasta`** (este passo). Confira:

```bash
cd ~/filogenomica
ls *.fasta
```

Devem aparecer os 8 arquivos.

## 4. Busca dos genes em todos os genomas

> **Antes de colar os comandos:** os comandos abaixo usam **todos os arquivos `.faa` da pasta**, então não é preciso digitar a lista de espécies. Confira só duas coisas com `ls *.faa`: (1) estão ali todos os proteômas escolhidos, e nenhum outro; (2) os nomes não têm espaços, por exemplo `Saccharomyces_cerevisiae_GCF_000146045.2.faa`. Esse nome (sem o `.faa`) será o nome da espécie na árvore.

### 4.1 Criar os bancos do BLAST

O BLAST por linha de comando precisa transformar cada genoma em um **banco de dados** antes de procurar nele. O comando `for` repete a mesma tarefa para cada arquivo `.faa` da pasta: a cada volta, a variável`$arquivo` é trocada pelo nome de um arquivo.

```bash
cd ~/filogenomica
for arquivo in *.faa; do
  makeblastdb -in $arquivo -dbtype prot -out $arquivo
done
```

### 4.2 Entendendo com um exemplo: uma proteína em um genoma

Antes de rodar tudo, faça o processo uma vez à mão. Nos comandos abaixo, troque `Kluyveromyces_lactis_GCF_000002515.2.faa` pelo nome de um dos seus arquivos `.faa` (veja com `ls *.faa`).

**Passo A — BLAST.** Procuramos a proteína RPO21 no genoma e guardamos só o nome do melhor resultado no arquivo `melhor_hit.txt`. As opções: `-evalue 1e-10` aceita só resultados muito significativos; `-max_target_seqs 1` e `-max_hsps 1` guardam só o melhor; `-outfmt "6 sseqid"` escreve apenas o nome da proteína encontrada.

```bash
cd ~/filogenomica
blastp -query RPO21.fasta \
  -db Kluyveromyces_lactis_GCF_000002515.2.faa \
  -evalue 1e-10 -max_target_seqs 1 -max_hsps 1 \
  -outfmt "6 sseqid" > melhor_hit.txt
cat melhor_hit.txt
```

Deve aparecer um código como `XP_...` ou `NP_...`: é a proteína desse genoma mais parecida com a RPO21.

**Obs**: O BLAST só escreve linhas na saída para queries que encontraram um hit. Portanto, se alguma sequência do arquivo `RPO21.fasta` não tiver alinhamento que satisfaça o limiar de e-value, ela **não aparecerá** no arquivo `melhor_hit.txt`

**Passo B — pegar a sequência.** O `seqkit grep` procura no genoma a proteína cujo nome está em `melhor_hit.txt (accession da proteina)` e mostra a sequência dela:

```bash
cd ~/filogenomica
seqkit grep -f melhor_hit.txt \
  Kluyveromyces_lactis_GCF_000002515.2.faa
```

**Passo C — trocar o nome e guardar.** O `|` envia a sequência para o `sed`, que troca a linha do nome (a que começa com `>`) pelo nome da espécie. O `>>` acrescenta o resultado no fim do arquivo da família, sem apagar o que já estava lá.

```bash
cd ~/filogenomica
seqkit grep -f melhor_hit.txt \
  Kluyveromyces_lactis_GCF_000002515.2.faa \
  | sed "s/^>.*/>Kluyveromyces_lactis_GCF_000002515.2/" \
  >> teste.fa
cat teste.fa
rm teste.fa
```

(O último comando apaga o arquivo de teste.)

### 4.3 Todas as proteínas em todos os genomas

Agora repetimos os passos A, B e C para todas as combinações, com dois loops `for`, um dentro do outro (executado dentro de um script):

Explicação dos comandos:

- `for proteina in *.fasta` passa por **cada proteína de referência** (todos os arquivos `.fasta` da pasta);
- `for genoma in *.faa` passa por **cada genoma** (todos os arquivos `.faa`).

Dois detalhes: `${genoma%.faa}` é o nome do arquivo do genoma sem o `.faa` (vira o nome da espécie), e `${proteina%.fasta}` é o nome da proteína sem o `.fasta` (vira o nome do arquivo da família, por exemplo `familia_RPO21.fa`). A primeira linha (`rm -f familia_*.fa`) apaga resultados antigos, para que o bloco possa ser repetido sem duplicar sequências.

Em vez de digitar esses comandos no terminal, vamos guardá-los em um **script**: um arquivo de texto com comandos, que o terminal executa de uma vez.

1. Abra o **Editor de Texto** do Ubuntu.
2. Copie e cole todo o conteúdo abaixo. A primeira linha (`#!/bin/bash`, chamada *shebang*) avisa que o arquivo deve ser lido pelo bash; as linhas que começam com `#` são comentários e não são executadas.

```bash
#!/bin/bash
# Busca cada proteina (.fasta) em cada genoma (.faa)
# e cria um arquivo familia_PROTEINA.fa por proteina.

rm -f familia_*.fa

for proteina in *.fasta; do
  for genoma in *.faa; do
    blastp -query $proteina -db $genoma \
      -evalue 1e-10 -max_target_seqs 1 -max_hsps 1 \
      -outfmt "6 sseqid" > melhor_hit.txt
    seqkit grep -f melhor_hit.txt $genoma \
      | sed "s/^>.*/>${genoma%.faa}/" \
      >> familia_${proteina%.fasta}.fa
  done
done

echo "Pronto!"
```

3. Salve o arquivo na pasta **`filogenomica`** com o nome **`busca_genes.sh`**.
4. No terminal, execute o script:

```bash
cd ~/filogenomica
bash busca_genes.sh
```

Quando aparecer `Pronto!`, a busca terminou. Para repetir a busca (por exemplo, depois de acrescentar proteínas), basta rodar `bash busca_genes.sh` de novo.

### 4.4 Conferir

```bash
cd ~/filogenomica
grep -c ">" familia_*.fa
cat familia_RPO21.fa
```

O primeiro comando conta as sequências de cada arquivo: deve ser igual ao número de espécies. O segundo mostra um dos arquivos: cada sequência começa com `>` e o nome da espécie com o código GCF.

Obs: os nomes das proteínas dentro do arquivo foram renomeadas com o nome da espécie para podermos fazer a combinação dos arquivos para a concatenação das sequências.

## 5. Alinhamento com o MAFFT

O MAFFT coloca os aminoácidos equivalentes de cada espécie na mesma coluna. O `for` passa por cada arquivo de família e cria um alinhamento com o mesmo nome, começando com `alinhado_`:

```bash
cd ~/filogenomica
for familia in familia_*.fa; do
  mafft --auto $familia > alinhado_$familia
done
```

Opcional: Para ver um alinhamento colorido, envie um arquivo `alinhado_....fa` ao [NCBI MSA Viewer](https://www.ncbi.nlm.nih.gov/projects/msaviewer/) (botão *Upload*).

## 6. Concatenação das sequências

Agora juntamos os 8 alinhamentos lado a lado em um único alinhamento grande, a **supermatriz**. O SeqKit junta as sequências que têm o mesmo nome; por isso, em todos os genes, cada sequência tem o mesmo nome do seu arquivo \`.faa\`.

```bash
cd ~/filogenomica
seqkit concat alinhado_*.fa > supermatriz.fa
seqkit stats supermatriz.fa
```

O `seqkit stats` mostra o número de espécies (`num_seqs`) e o comprimento do alinhamento (`max_len`).

## 7. Cálculo da árvore com o IQ-TREE

```bash
cd ~/filogenomica
iqtree -s supermatriz.fa -m MFP -B 1000
```

O que cada opção faz:

- `-s supermatriz.fa`: o alinhamento de entrada.
- `-m MFP`: o IQ-TREE testa vários modelos de evolução de proteínas e escolhe sozinho o melhor (ModelFinder).
- `-B 1000`: calcula o suporte de cada ramo com 1000 réplicas de *ultrafast bootstrap *
- Alternativa: você pode rodar o iqtree indicando o grupo externo.&#32;

  ```bash
  iqtree -s supermatriz.fa -m MFP -B 1000 \
    -o "nome exato do grupo externo na supermatriz"
  ```

**Qual modelo foi escolhido?**

```bash
cd ~/filogenomica
grep "Best-fit model" supermatriz.fa.iqtree
```

Anote o modelo (por exemplo, `LG+F+I+G4`). O comando "grep" buscou a linha que informa o modelo escolhido automaticamente pelo iqtree. O arquivo também pode ser lido com um editor de texto.

**Ver a árvore:** abra [itol.embl.de](https://itol.embl.de), clique em **Upload** e envie o arquivo `supermatriz.fa.treefile`, que está na pasta `filogenomica`. No painel do iTOL, ative a exibição dos valores de bootstrap.

## Se algo der errado

| Mensagem | O que fazer |
| --- | --- |
| `command not found` | Rode `micromamba activate filogenomica`. |
| `No such file or directory` | Rode `cd ~/filogenomica` e `ls`. Confira se os nomes dos arquivos estão exatamente iguais aos do roteiro (maiúsculas, `.faa` e `.fasta`). |
| Um `familia_*.fa` tem menos sequências que o esperado | O BLAST não achou o gene em alguma espécie. Mova o arquivo .fasta dessa proteína para fora da pasta (ou apague-o) e repita a partir do passo 4.3. |
| `mkdir: cannot create directory: File exists` | A pasta já existe. Pode ignorar. |
