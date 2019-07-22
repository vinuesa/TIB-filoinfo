# <b>TIBS-filoinfo</b>
### <b>Talleres Internacionales de Bioinformática - Centro de Ciencias Genómicas, UNAM, Cuernavaca, México</b>

Este repositorio contiene el material para el [Taller 3 - An&aacute;lisis comparativo de genomas microbianos: Pangen&oacute;mica y filoinform&aacute;tica](http://congresos.nnb.unam.mx/TIB2019/t3-analisis-comparativo-de-genomas-microbianos-pangenomica-y-filoinformatica/) de los [Talleres Internacionales de Bioinform&aacute;tica - TIB2019](http://congresos.nnb.unam.mx/TIB2019), celebrados en el [Centro de Ciencias Genómicas](http://www.ccg.unam.mx) de la [Universidad Nacional Aut&oacute;noma de M&eacute;xico](http://www.ccg.unam.mx), del 29 de julio al 2 de agosto de 2019.

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
Es necesario que el alumno traiga su computadora personal, de preferencia con Linux (o MacOS X) como sistema operativo. Si usan Windows, deberán tener instalado [MobaXterm](https://mobaxterm.mobatek.net/download-home-edition.html) (para Ms Windows), o alguna otra herramienta que permita establecer una conexión SSH al servidor que corre los programas a usar en el curso.

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
Si tienes instalado [git](https://git-scm.com/) en tu computadora, puedes clonar el repositorio con el comando:

   <code>git clone https://github.com/vinuesa/TIBS-filoinfo.git</code>

En [ubuntu](https://www.ubuntu.com/) es muy f&aacute;cil instalar git: 

  <code>sudo apt install git</code>

***

## Sesiones y material asociado

### Horario y lugar de impartici&oacute;n de las sesiones
Las clases se imparten del 29 de Julio al 2 de Agosto en el aula 3 de la LCG-UNAM, Cuernavaca, Morelos
de 9 a 17:30 hrs, seg&uacute;n el [programa de los TIB2019](http://congresos.nnb.unam.mx/TIB2019/programa/)

#### Sesiones
- Introducción a Linux (teoría y práctica)
- Conceptos básicos de biología evolutiva y filogenética
- Búsqueda de homólogos usando BLAST desde la línea de comandos (prácticas)
- Alineamientos múltiples (prácticas)
- Introducción a los métodos filogenéticos, árboles de genes y de árboles de especies
- Modelos de sustitución y máxima verosimilitud (teoría)
- Ajuste de modelos e inferencia de filogenias de máxima verosimilitud (prácticas)
- Delimitación de especies bacterianas usando métodos evolutivos y datos multilocus
- Inferencia bayesiana de filogenias (teoría y práctica)
- Pangenómica y evolución microbiana (Seminario de investigación)
- Cómputo de familias de genes homólogos con datos genómicos (teoría)
- Análisis pangenómico usando GET_HOMOLOGUES (prácticas)
- Estrategias para la estima de filogenias genómicas (teoría)
- Estima de filogenias genómicas con GET_PHYLOMARKERS (prácticas)

En construcci&oacute;n ...



