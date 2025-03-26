pro iris_test2
;这个程序是用来看如何把IDL的数据存储下来再让PYTHON读出来
a=[[[2,2],[2.5,2.5],[2,2],[2,2]],[[3,3],[3,3],[3,3],[3,3]],[[4,4],[4,4],[4,4],[4,4]]]
help,a
a=float(a)
print,a
filename = 'matrix_test.dat'
matrix=a
; 打开文件以二进制写入
OPENW, unit, filename, /BINARY, /GET_LUN  ; 打开文件以二进制写入

; 将矩阵数据写入文件
; 注意：这里我们直接将整个矩阵作为一个大的数据块写入文件
; 在Python中读取这个文件时，你需要知道数据的维度和类型来正确地解析它
WRITEU, unit, matrix

; 关闭文件
FREE_LUN, unit
end