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

#### <b>Sesión 1: Introducción a Linux (teoría y práctica)</b>
- [presentaci&oacute;n - PDF: Primer contacto con un sistema GNU/Linux](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion1_intro2linux/Intro_biocomputo_Linux_pt1.pdf)
- Pr&aacute;ctica. Navegación del sistema, uso de comandos básicos y ejercicio de parseo de archivo FASTA
      + [pr&aacute;ctica - html](https://vinuesa.github.io/TIB-filoinfo/sesion1_intro2linux/) 
      + [pr&aacute;ctica - pdf](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion1_intro2linux/working_with_linux_commands.pdf)
  
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


#### <b>Sesión 2: Conceptos básicos de biología evolutiva y filogenética</b>
- [presentaci&oacute;n - PDF: conceptos básicos de filogenética y evolución](https://github.com/vinuesa/TIB-filoinfo/tree/master/docs/sesion2_conceptos_basicos/sesion2_conceptos_basicos_filogenetica_y_evol.pdf)

#### <b>Sesión 3: Búsqueda de homólogos usando BLAST desde la línea de comandos (teoría y prácticas)</b>
#### <b>Sesión 4: Alineamientos múltiples (teoría y prácticas)</b>
#### <b>Sesión 5: Introducción a los métodos filogenéticos, árboles de genes y de árboles de especies</b>
#### <b>Sesión 6: Modelos de sustitución y máxima verosimilitud (teoría)</b>
#### <b>Ajuste de modelos e inferencia de filogenias de máxima verosimilitud (prácticas)</b>
#### <b>Delimitación de especies bacterianas usando métodos evolutivos y datos multilocus</b>
#### <b>Inferencia bayesiana de filogenias (teoría y práctica)</b>
#### <b>Pangenómica y evolución microbiana (Seminario de investigación)</b>
#### <b>Cómputo de familias de genes homólogos con datos genómicos (teoría)</b>
#### <b>Análisis pangenómico usando GET_HOMOLOGUES (prácticas)</b>
#### <b>Estrategias para la estima de filogenias genómicas (teoría)</b>
#### <b>Estima de filogenias genómicas con GET_PHYLOMARKERS (prácticas)</b>

En construcci&oacute;n ...



