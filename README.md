# Computer Methods for Mathematical Computations

by Nikolai Shokhirev

### 1. Introduction

Computer Methods for Mathematical Computations by George Forsythe, Michael Malcolm and Cleve Moler (FMM) is one of the great classic textbooks on numerical methods. First I read it in Russian translation [1] and later bought the original [2].

I used the book both as a reference to numerical concepts and the source of working efficient algorithms. During my work as a researcher and scientific programmer, I translatedand the original FORTRAN code to several languages.

I am going to share these translations on GitHub: FMM ( Direct link: https://github.com/dr-nikolai/FMM ). All updates will be documented at FMM Wiki..

### 2. Fortran-90 translation

The original FMM algorithms are perfectly usable and available at NetLib. However, Ralph Carmichael made translation to a modern version of FOFTRAN: [3]. He also added several testing routines.

### 3. Pascal translation

I translated the FMM routines to Object Pascal (Delphi). They are included into my numerical library PasMatLib [4].

These algorithm, along with some dependencies and DUnit tests are collected in fmm.pas [5]. They can be used independently from PasMatLib.



### References

 1. &#1060;&#1086;&#1088;&#1089;&#1072;&#1081;&#1090; &#1044;&#1078;., &#1052;&#1072;&#1083;&#1100;&#1082;&#1086;&#1083;&#1100;&#1084; &#1052;.,
  &#1052;&#1086;&#1091;&#1083;&#1077;&#1088; &#1050;. &#1052;&#1072;&#1096;&#1080;&#1085;&#1085;&#1099;&#1077; &#1084;&#1077;&#1090;&#1086;&#1076;&#1099; &#1084;&#1072;&#1090;&#1077;&#1084;&#1072;&#1090;&#1080;&#1095;&#1077;&#1089;&#1082;&#1080;&#1093; 
  &#1074;&#1099;&#1095;&#1080;&#1089;&#1083;&#1077;&#1085;&#1080;&#1081;. &#1052;.: &#1052;&#1080;&#1088;, 1980.
 2. Forsythe, George E.; Malcolm, Michael A.; Moler, Cleve B.: Computer Methods for Mathematical Computations, Prentice-Hall, 1977.
 3. [FMM Fortran-90](http://www.pdas.com/fmmdownload.html) - Computer Methods for Mathematical Computations. 
 4. [PasMatLib](http://www.shokhirev.com/nikolai/abc/sciprog/DTutorial.html)- Scientific programming with Delphi. 
 5. [fmm.pas](https://github.com/dr-nikolai/FMM/tree/master/fmm.pas) - Delphi/Pascal implementation (Delphi 5+, Lazarus) by Nikolai Shokhirev

