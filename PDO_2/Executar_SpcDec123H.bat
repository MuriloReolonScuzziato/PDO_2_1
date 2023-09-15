REM ----------------------
REM SpcDec3 deterministico
REM ----------------------
copy "D:\Murilo\PDO\det\*.osi" "D:\Murilo\PDO\A"
cd D:\Murilo\PDO\A
del report.txt
rename D:\Murilo\PDO\Dados\Sistema_46bar\A\Cond_Inic_H_m.txt Cond_Inic_H.txt
for /l %%x in (1, 1, 30) do (
   PDO_SpcDec3H 2 1 0 10 1 1 1 A/ 0 %%x 1 1 1
)
move log_* D:\Murilo\PDO\Resultados\Sistema_46bar\A\Hrst\SpcDec3
rename report.txt reportAm.txt
move reportAm.txt D:\Murilo\PDO\Resultados\Sistema_46bar\A\Hrst\SpcDec3

copy "D:\Murilo\PDO\det\*.osi" "D:\Murilo\PDO\B"
cd D:\Murilo\PDO\B
del report.txt
rename D:\Murilo\PDO\Dados\Sistema_46bar\B\Cond_Inic_H_m.txt Cond_Inic_H.txt
for /l %%x in (1, 1, 30) do (
   PDO_SpcDec3H 2 1 0 10 1 1 1 B/ 0 %%x 1 1 1
)
move log_* D:\Murilo\PDO\Resultados\Sistema_46bar\B\Hrst\SpcDec3
rename report.txt reportBm.txt
move reportBm.txt D:\Murilo\PDO\Resultados\Sistema_46bar\B\Hrst\SpcDec3

copy "D:\Murilo\PDO\det\*.osi" "D:\Murilo\PDO\C"
cd D:\Murilo\PDO\C
del report.txt
rename D:\Murilo\PDO\Dados\Sistema_46bar\C\Cond_Inic_H_m.txt Cond_Inic_H.txt
for /l %%x in (1, 1, 30) do (
   PDO_SpcDec3H 2 1 0 10 1 1 1 C/ 0 %%x 1 1 1
)
move log_* D:\Murilo\PDO\Resultados\Sistema_46bar\C\Hrst\SpcDec3
rename report.txt reportCm.txt
move reportCm.txt D:\Murilo\PDO\Resultados\Sistema_46bar\C\Hrst\SpcDec3

copy "D:\Murilo\PDO\det\*.osi" "D:\Murilo\PDO\D"
cd D:\Murilo\PDO\D
del report.txt
rename D:\Murilo\PDO\Dados\Sistema_46bar\D\Cond_Inic_H_m.txt Cond_Inic_H.txt
for /l %%x in (1, 1, 30) do (
   PDO_SpcDec3H 2 1 0 10 1 1 1 D/ 0 %%x 1 1 1
)
move log_* D:\Murilo\PDO\Resultados\Sistema_46bar\D\Hrst\SpcDec3
rename report.txt reportDm.txt
move reportDm.txt D:\Murilo\PDO\Resultados\Sistema_46bar\D\Hrst\SpcDec3

copy "D:\Murilo\PDO\det\*.osi" "D:\Murilo\PDO\E"
cd D:\Murilo\PDO\E
del report.txt
rename D:\Murilo\PDO\Dados\Sistema_46bar\E\Cond_Inic_H_m.txt Cond_Inic_H.txt
for /l %%x in (1, 1, 30) do (
   PDO_SpcDec3H 2 1 0 10 1 1 1 E/ 0 %%x 1 1 1
)
move log_* D:\Murilo\PDO\Resultados\Sistema_46bar\E\Hrst\SpcDec3
rename report.txt reportEm.txt
move reportEm.txt D:\Murilo\PDO\Resultados\Sistema_46bar\E\Hrst\SpcDec3


REM ----------------------
REM SpcDec2 deterministico
REM ----------------------
cd D:\Murilo\PDO\A
for /l %%x in (1, 1, 30) do (
   PDO_SpcDec2H 2 1 0 10 1 1 1 A/ 0 %%x 1 1 1
)
move log_* D:\Murilo\PDO\Resultados\Sistema_46bar\A\Hrst\SpcDec2
rename report.txt reportAm.txt
move reportAm.txt D:\Murilo\PDO\Resultados\Sistema_46bar\A\Hrst\SpcDec2

cd D:\Murilo\PDO\B
for /l %%x in (1, 1, 30) do (
   PDO_SpcDec2H 2 1 0 10 1 1 1 B/ 0 %%x 1 1 1
)
move log_* D:\Murilo\PDO\Resultados\Sistema_46bar\B\Hrst\SpcDec2
rename report.txt reportBm.txt
move reportBm.txt D:\Murilo\PDO\Resultados\Sistema_46bar\B\Hrst\SpcDec2

cd D:\Murilo\PDO\C
for /l %%x in (1, 1, 30) do (
   PDO_SpcDec2H 2 1 0 10 1 1 1 C/ 0 %%x 1 1 1
)
move log_* D:\Murilo\PDO\Resultados\Sistema_46bar\C\Hrst\SpcDec2
rename report.txt reportCm.txt
move reportCm.txt D:\Murilo\PDO\Resultados\Sistema_46bar\C\Hrst\SpcDec2

cd D:\Murilo\PDO\D
for /l %%x in (1, 1, 30) do (
   PDO_SpcDec2H 2 1 0 10 1 1 1 D/ 0 %%x 1 1 1
)
move log_* D:\Murilo\PDO\Resultados\Sistema_46bar\D\Hrst\SpcDec2
rename report.txt reportDm.txt
move reportDm.txt D:\Murilo\PDO\Resultados\Sistema_46bar\D\Hrst\SpcDec2

cd D:\Murilo\PDO\A
for /l %%x in (1, 1, 30) do (
   PDO_SpcDec2H 2 1 0 10 1 1 1 E/ 0 %%x 1 1 1
)
move log_* D:\Murilo\PDO\Resultados\Sistema_46bar\E\Hrst\SpcDec2
rename report.txt reportEm.txt
move reportEm.txt D:\Murilo\PDO\Resultados\Sistema_46bar\E\Hrst\SpcDec2


REM ----------------------
REM SpcDec1 deterministico
REM ----------------------
cd D:\Murilo\PDO\A
for /l %%x in (1, 1, 30) do (
   PDO_SpcDec1H 2 1 0 10 1 1 1 A/ 0 %%x 1 1 1
)
del *.osi
move log_* D:\Murilo\PDO\Resultados\Sistema_46bar\A\Hrst\SpcDec1
rename report.txt reportAm.txt
move reportAm.txt D:\Murilo\PDO\Resultados\Sistema_46bar\A\Hrst\SpcDec1
rename D:\Murilo\PDO\Dados\Sistema_46bar\A\Cond_Inic_H.txt Cond_Inic_H_m.txt

cd D:\Murilo\PDO\B
for /l %%x in (1, 1, 30) do (
   PDO_SpcDec1H 2 1 0 10 1 1 1 B/ 0 %%x 1 1 1
)
del *.osi
move log_* D:\Murilo\PDO\Resultados\Sistema_46bar\B\Hrst\SpcDec1
rename report.txt reportBm.txt
move reportBm.txt D:\Murilo\PDO\Resultados\Sistema_46bar\B\Hrst\SpcDec1
rename D:\Murilo\PDO\Dados\Sistema_46bar\B\Cond_Inic_H.txt Cond_Inic_H_m.txt

cd D:\Murilo\PDO\C
for /l %%x in (1, 1, 30) do (
   PDO_SpcDec1H 2 1 0 10 1 1 1 C/ 0 %%x 1 1 1
)
del *.osi
move log_* D:\Murilo\PDO\Resultados\Sistema_46bar\C\Hrst\SpcDec1
rename report.txt reportCm.txt
move reportCm.txt D:\Murilo\PDO\Resultados\Sistema_46bar\C\Hrst\SpcDec1
rename D:\Murilo\PDO\Dados\Sistema_46bar\C\Cond_Inic_H.txt Cond_Inic_H_m.txt

cd D:\Murilo\PDO\D
for /l %%x in (1, 1, 30) do (
   PDO_SpcDec1H 2 1 0 10 1 1 1 D/ 0 %%x 1 1 1
)
del *.osi
move log_* D:\Murilo\PDO\Resultados\Sistema_46bar\D\Hrst\SpcDec1
rename report.txt reportDm.txt
move reportDm.txt D:\Murilo\PDO\Resultados\Sistema_46bar\D\Hrst\SpcDec1
rename D:\Murilo\PDO\Dados\Sistema_46bar\D\Cond_Inic_H.txt Cond_Inic_H_m.txt

cd D:\Murilo\PDO\E
for /l %%x in (1, 1, 30) do (
   PDO_SpcDec1H 2 1 0 10 1 1 1 E/ 0 %%x 1 1 1
)
del *.osi
move log_* D:\Murilo\PDO\Resultados\Sistema_46bar\E\Hrst\SpcDec1
rename report.txt reportEm.txt
move reportEm.txt D:\Murilo\PDO\Resultados\Sistema_46bar\E\Hrst\SpcDec1
rename D:\Murilo\PDO\Dados\Sistema_46bar\E\Cond_Inic_H.txt Cond_Inic_H_m.txt