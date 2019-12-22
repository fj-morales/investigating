# Investigating the case of weak baselines in Ad-Hoc retrieval and Question Answering repo

Accompanying repository for my MSc thesis: "Investigating the case of weak baselines in Ad-Hoc retrieval and Question Answering"

For the code accompanying paper A, "Deep Relevance Ranking Using Enhanced Document-Query Interactions". [[PDF](http://nlp.cs.aueb.gr/pubs/emnlp2018.pdf)]], click [[here](https://github.com/fj-morales/deep-relevance-ranking)]

For the code accompanying paper B, "TVQA: Localized, Compositional Video Question Answering". [[PDF](https://www.aclweb.org/anthology/D18-1167.pdf)], ckick [[here](https://github.com/fj-morales/TVQA)]

## Prerequisites

The execution of these repositories depend on the next prerequisites:

**Oracle Java 11**

Make sure you have a running insatllation of Oracle Java 11.

`sudo apt-get install oracle-java11-installer-local`

**RankLib 2.12**

Get a copy of the .jar file from [[Ranklib 2.12](https://sourceforge.net/projects/lemur/files/lemur/RankLib-2.12/RankLib-2.12.jar/download)] and keep it in the directory `./ranklib/RankLib-2.12.jar` at the same parent level as the papers' repos.

**trec_eval**

Compile the executable of [[trec_eval](https://github.com/usnistgov/trec_eval)] and keep it in the directory `./trec_eval/trec_eval` at the same parent level as the papers' repos.

**HpBandSter**

Install HpBandSter for HPO (only compatible with Python 3!):

`pip install hpbandster`
