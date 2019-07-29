# <img src="docs/pics/NNB-TIB-Logo.png" alt="TIB2019" width="115" height="115" align="middle"> <b>TIB2019-T3</b>

<img src="docs/pics/TIB2019_banner.png" alt="TIB2019-screenshot">

## <b>Talleres Internacionales de Bioinformática - Centro de Ciencias Genómicas, UNAM, Cuernavaca, México</b>

### Sobre el repositorio
Este repositorio contiene el material para el [Taller 3 - An&aacute;lisis comparativo de genomas microbianos: Pangen&oacute;mica y filoinform&aacute;tica](http://congresos.nnb.unam.mx/TIB2019/t3-analisis-comparativo-de-genomas-microbianos-pangenomica-y-filoinformatica/) de los [Talleres Internacionales de Bioinform&aacute;tica - TIB2019](http://congresos.nnb.unam.mx/TIB2019), a celebrarse en el [Centro de Ciencias Genómicas](http://www.ccg.unam.mx) de la [Universidad Nacional Aut&oacute;noma de M&eacute;xico](http://www.ccg.unam.mx), del 29 de julio al 2 de agosto de 2019.

#### Otros repositorios asociados a ediciones anteriores de los TIB
- [T2: Análisis exploratorio y estadístico de datos biológicos usando R](https://github.com/vinuesa/curso_Rstats), edición [TIB2018](http://congresos.nnb.unam.mx/TIB2018/)

***
 
## Presentaci&oacute;n
Este taller de nivel básico-intermedio te proporcionará una sólida base de conocimientos teóricos y prácticos sobre aspectos fundamentales de biocómputo para inferencia filogenética, evolución molecular y genómica microbiana, con énfasis en análisis pangenómicos y filogenómicos.

### Descripción
En el taller (40 hrs) tendremos sesiones teóricas y prácticas que cubrirán un amplio espectro de aspectos básicos del tópico como:

- escrutinio de bases de datos mediante BLAST
- determinación e interpretación de homología
- alineamiento de múltiples secuencias y conversión de formatos 
- inferencia filogenética
- análisis pangenómico y filogenómico de genomas microbianos

Se darán presentaciones detalladas del uso de programas clave (todos de “open source”) para estos análisis, usando datos tomados de las bases de datos. También se presentará el uso de algunos scripts de Bash y Perl muy sencillos, con el objetivo de aprender los aspectos básicos de estos lenguajes para el análisis de datos genómicos.

Al final del curso tendrán una amplia visión sobre el espectro de posibilidades que brindan la filogenética y la evolución molecular en distintos tipos de estudios biológicos y genómicos, que les servirán como herramientas conceptuales y metodológicas de gran utilidad en su carrera como estudiantes o profesionales.

### Requisitos
#### Conocimientos previos
Es recomendable tener conocimientos básicos de Unix/Linux a nivel básico, ya que todas las demostraciones de software se harán en este sistema operativo.

#### Requisitos técnicos
Es necesario que el alumno traiga su computadora personal, de preferencia con Linux (o MacOS X) como sistema operativo. 

<b>Si usan Windows, deberán tener instalado [MobaXterm](https://mobaxterm.mobatek.net/download-home-edition.html) (para Ms Windows) antes de llegar al taller!</b>. 

Aquí tienen [instrucciones para la instalación de MobaXterm en Windows](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/Instalación_de_mobaXterm_en_Windows.pdf)

### Sobre el profesor
Hola, me llamo [Pablo Vinuesa](http://www.ccg.unam.mx/~vinuesa/). Soy investigador titular del 
[Centro de Ciencias Gen&oacute;micas](http://www.ccg.unam.mx) de la 
[Universidad Nacional Aut&oacute;noma de M&eacute;xico - UNAM](http://www.unam.mx/).

Mis [l&iacute;neas de investigaci&oacute;n](http://www.ccg.unam.mx/~vinuesa/research.html) 
integran la gen&oacute;mica y la bioinform&aacute;tica con la biolog&iacute;a y gen&eacute;tica molecular para entender 
la evoluci&oacute;n y emergencia de pat&oacute;genos oportunistas a partir de microbios ambientales.

### Sobre el material did&aacute;ctico
A trav&eacute;s de estas p&aacute;ginas se distribuyen los apuntes, ejercicios y datos que se usar&aacute;n en el [Taller 3 - An&aacute;lisis comparativo de genomas microbianos: Pangen&oacute;mica y filoinform&aacute;tica](http://congresos.nnb.unam.mx/TIB2019/t3-analisis-comparativo-de-genomas-microbianos-pangenomica-y-filoinformatica/).
Para tu convenienca, se distribuye en formatos pdf y html.

Puedes ver en mi sitio Web el [listados de cursos](http://www.ccg.unam.mx/~vinuesa/cursos.html) y materiales asociados, que pongo libremente disponible para la comunidad.

### Licencia y términos de uso
El material del [T3, TIB-filoinfo](http://congresos.nnb.unam.mx/TIB2019/t3-analisis-comparativo-de-genomas-microbianos-pangenomica-y-filoinformatica/) lo distribuyo p&uacute;blicamente a trav&eacute;s de este repositorio GitHub bajo la [**Licencia No Comercial Creative Commons 4.0**](https://creativecommons.org/licenses/by-nc/4.0/) 

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0</a>


### Clonaci&oacute;n del repositorio
- Si tienes instalado [git](https://git-scm.com/) en tu computadora, puedes clonar el repositorio con el comando:

   <code>git clone https://github.com/vinuesa/TIBS-filoinfo.git</code>

- Para actualizar el repositorio, ejecuta este comando desde dentro del directorio TIBS-filoinfo
  
   <code>git pull https://github.com/vinuesa/TIBS-filoinfo.git</code>

En [ubuntu](https://www.ubuntu.com/) y [MobaXterm](https://mobaxterm.mobatek.net/download-home-edition.html) es muy f&aacute;cil instalar git: 

  <code>sudo apt install git</code>

Vean además las [instrucciones para la instalación de MobaXterm en Windows](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/Instalación_de_mobaXterm_en_Windows.pdf), 
que indican cómo instalar el <i>Git plugin</i> de MobaXterm. 

***

## <b>Sesiones y material asociado</b>
### Horario y lugar de impartici&oacute;n de las sesiones
Las clases se imparten del 29 de Julio al 2 de Agosto en el auditorio Guillermo Soberón del [CCG-UNAM]((http://www.ccg.unam.mx/), Cuernavaca, Morelos
de 9 a 17:30 hrs, seg&uacute;n el [programa de los TIB2019](http://congresos.nnb.unam.mx/TIB2019/programa/)


### Lunes 29 de Julio

#### <b>Sesión 1: Introducción a Linux (teoría y práctica)</b>
- [presentaci&oacute;n - PDF: Primer contacto con un sistema GNU/Linux](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion1_intro2linux/Intro_biocomputo_Linux_pt1.pdf)
- Pr&aacute;ctica. Navegación del sistema, uso de comandos básicos y ejercicio de parseo de archivo FASTA
  - [pr&aacute;ctica - html](https://vinuesa.github.io/TIB-filoinfo/sesion1_intro2linux/) 
  - [pr&aacute;ctica - pdf](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion1_intro2linux/working_with_linux_commands.pdf)
  
<!--
+ Pr&aacute;ctica 2. Descarga de secuencias en formato FASTA de GenBank usando el sistema ENTREZ y parseo de los archivos usando herrramientas de filtrado
    - [pr&aacute;ctica2 - html](https://vinuesa.github.io/OMICAS_UAEM/practica2_parseo_fastas/)
    - [pr&aacute;ctica2 - pdf](https://github.com/vinuesa/OMICAS_UAEM/tree/master/docs/practica2_parseo_fastas/ejercicio_parseo_fastas_ENTREZ.pdf)
-->

- Lecturas recomendadas:
  - Atma Ivancevic. The ten commandments for learning how to code. [Carrer Column, Nature, 20 Feb. 2019](https://www.nature.com/articles/d41586-019-00653-5)
  - Velez Rueda AJ, Benítez GI, Marchetti J, Hasenahuer MA, Fornasari MS, Palopoli N, Parisi G. Bioinformatics calls the school: Use of smartphones to introduce
Python for bioinformatics in high schools. [PLoS Comput Biol. 2019 Feb 14;15(2):e1006473.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006473)
  - Dudley JT, Butte AJ. A quick guide for developing effective bioinformatics programming skills. [PLoS Comput Biol. 2009 Dec;5(12):e1000589](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000589)
  - The Linux Command Line - a complete introduction. [William E. Shotts, Jr. No Starch Press](http://linuxcommand.org/lc3_learning_the_shell.php#contents)
  - Bioinformatics Data Skills: Reproducible and Robust Research with Open Source Tools. [Vince Buffalo. O'Reilly Media 2014](http://freecomputerbooks.com/Bioinformatics-Data-Skills.html)

#### <b>Sesión 2: Conceptos básicos de biología evolutiva, filogenética y (pan)genómica microbiana</b>
- [presentaci&oacute;n - PDF: conceptos básicos de filogenética y evolución](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion2_conceptos_basicos/sesion2_conceptos_basicos_filogenetica_y_evol.pdf)

- Lecturas recomendadas
  - Fitch WM. Homology a personal view on some of the problems. [Trends Genet. 2000 May;16(5):227-31](https://www.ncbi.nlm.nih.gov/pubmed/10782117)
  - Koonin EV. Orthologs, paralogs, and evolutionary genomics. [Annu Rev Genet. 2005;39:309-38](https://www.ncbi.nlm.nih.gov/pubmed/16285863)
  - Glover N, Dessimoz C, Ebersberger I, Forslund SK, Gabaldón T, Huerta-Cepas J, Martin MJ et al. Quest for Orthologs Consortium. Advances and Applications in the Quest for Orthologs. [Mol Biol Evol. 2019 Jun 26. pii: msz150. doi: 10.1093/molbev/msz150.](https://www.ncbi.nlm.nih.gov/pubmed/31241141)
  - Vernikos G, Medini D, Riley DR, Tettelin H. Ten years of pan-genome analyses. [Curr Opin Microbiol. 2015 Feb;23:148-54](https://www.ncbi.nlm.nih.gov/pubmed/25483351)
  - McInerney JO, McNally A, O'Connell MJ. Why prokaryotes have pangenomes. [Nat Microbiol. 2017 Mar 28;2:17040](https://www.ncbi.nlm.nih.gov/pubmed/28350002)
  - Sela I, Wolf YI, Koonin EV. Theory of prokaryotic genome evolution. [Proc Natl Acad Sci U S A. 2016 Oct 11;113(41):11399-11407](https://www.ncbi.nlm.nih.gov/pubmed/27702904)
  - Land M, Hauser L, Jun SR, Nookaew I, Leuze MR, Ahn TH, Karpinets T, Lund O, Kora G, Wassenaar T, Poudel S, Ussery DW. Insights from 20 years of bacterial genome sequencing. [Funct Integr Genomics. 2015 Mar;15(2):141-61](https://www.ncbi.nlm.nih.gov/pubmed/25722247)

***

### Martes 30 de Julio

#### <b>Sesión 3: Búsqueda de homólogos usando BLAST desde la línea de comandos (teoría y prácticas)</b>
- [presentación - PDF](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion3_BLAST/Tema3_BLAST_OVERVIEW.pdf)
- práctica
  - [comandos, txt](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion3_BLAST/data/runing_and_parsing_BLAST_from_the_cmd_line.commands)
  - [16S_4blastN.tgz](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion3_BLAST/data/16S_4blastN.tgz)
  - [gene_discovery_and_annotation_using_blastx.tgz](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion3_BLAST/data/gene_discovery_and_annotation_using_blastx.tgz)
  - [split_fasta.pl](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/bin/split_fasta.pl)
  - [blast-imager.pl](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/bin/blast-imager.pl)
- Lecturas recomendadas
  - Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. [BLAST+: architecture and applications. BMC Bioinformatics. 2009 Dec 15;10:421](https://www.ncbi.nlm.nih.gov/pubmed/20003500)
  - Hu G, Kurgan L. Sequence Similarity Searching. [Curr Protoc Protein Sci. 2019 Feb;95(1):e71. doi: 10.1002/cpps.71](https://www.ncbi.nlm.nih.gov/pubmed/30102464)

#### <b>Sesión 4: Alineamientos múltiples (teoría y prácticas)</b>
- [presentación - PDF](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion4_alineamientos/Tema4_alineamientos_multiples.pdf)
- práctica
  - [comandos, txt](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion4_alineamientos/practicas_aln_multiples_clustal.cmds)
  - [sequences, tgz](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion4_alineamientos/sequences_for_alingment.tgz)
  - [align_seqs_with_clustal_or_muscle.sh](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/bin/align_seqs_with_clustal_or_muscle.sh)
  - [convert_alnFormats_using_clustalw.sh](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/bin/convert_alnFormats_using_clustalw.sh)
  - [convert_aln_format_batch_bp.pl](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/bin/convert_aln_format_batch_bp.pl)
  - [translate_fastas.pl](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/bin/translate_fastas.pl)
  - [prot2cdnAlns.pl](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/bin/prot2cdnAlns.pl)
- Lecturas recomendadas
  - Simossis V, Kleinjung J, Heringa J. An overview of multiple sequence alignment. [Curr Protoc Bioinformatics. 2003 Nov;Chapter 3:Unit 3.7](https://www.ncbi.nlm.nih.gov/pubmed/18428699)
  - Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG. Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. [Mol Syst Biol. 2011 Oct 11;7:539](https://www.ncbi.nlm.nih.gov/pubmed/21988835)
  - Sievers F, Higgins DG. Clustal Omega for making accurate alignments of many protein sequences. [Protein Sci. 2018 Jan;27(1):135-145](https://www.ncbi.nlm.nih.gov/pubmed/28884485)

#### <b>Sesión 5: Introducción a los métodos filogenéticos, modelos de sustitución y algoritmos de búsqueda de árboles</b>
- [presentación - PDF](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion5_metodos_filogeneticos_y_modelos/Tema5_metodos_filogeneticos_y_modelos.pdf)
- Lecturas recomendadas
  - Yang Z, Rannala B. Molecular phylogenetics: principles and practice. [Nat Rev Genet. 2012 Mar 28;13(5):303-14](https://www.ncbi.nlm.nih.gov/pubmed/22456349)

### Miércoles 31 de Julio

***

#### <b>Sesión 6: Selección de modelos e inferencia de filogenias bajo máxima verosimilitud (teoría y práctica)</b>
- [presentación - PDF](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion6_maxima_verosimilitud/Tema6_maxima_verosimilitud_y_seleccion_de_models.pdf)
- práctica
  - [tutorial phyml, comandos - html](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion6_maxima_verosimilitud/)
  - [tutorial phyml (secuencias), tgz](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion6_maxima_verosimilitud/phyml_tutorial_data.tgz)
  - [tutorial jmodeltest (comandos y secuencias), tgz](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion6_maxima_verosimilitud/jmodeltest_tutorial.tgz)
- Lecturas recomendadas
  - Lefort V, Longueville JE, Gascuel O. SMS: Smart Model Selection in PhyML. [Mol Biol Evol. 2017 Sep 1;34(9):2422-2424](https://www.ncbi.nlm.nih.gov/pubmed/28472384)
  - Criscuolo A. morePhyML: improving the phylogenetic tree space exploration with PhyML 3. Mol [Phylogenet Evol. 2011 Dec;61(3):944-8](https://www.ncbi.nlm.nih.gov/pubmed/21925283)
  - Guindon S, Dufayard JF, Lefort V, Anisimova M, Hordijk W, Gascuel O. New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. [Syst Biol. 2010 May;59(3):307-21](https://www.ncbi.nlm.nih.gov/pubmed/20525638)


#### <b>Inferencia bayesiana de filogenias (teoría y práctica)</b>
- Lecturas recomendadas
  - Nascimento FF, Reis MD, Yang Z. A biologist's guide to Bayesian phylogenetic analysis. [Nat Ecol Evol. 2017 Oct;1(10):1446-1454](https://www.ncbi.nlm.nih.gov/pubmed/28983516)

***


### Jueves 1 de Agosto
#### <b>Pangenómica y evolución microbiana (Seminario de investigación)</b>
#### <b>Cómputo de familias de genes homólogos con datos genómicos (teoría)</b>
#### <b>Análisis pangenómico usando GET_HOMOLOGUES (prácticas)</b>
- Lecturas recomendadas
  -
***

### Viernes 2 de Agosto
#### <b>Estrategias para la estima de filogenias genómicas (teoría)</b>
#### <b>Estima de filogenias genómicas con GET_PHYLOMARKERS (prácticas)</b>
- Lecturas recomendadas
  -
En construcci&oacute;n ...



