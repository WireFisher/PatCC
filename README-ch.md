# PatCC1使用手册

PatCC1 (**Pa**rallel **T**riangulation Algorithm with **C**ommonality and Parallel **C**onsistency, Version 1)

## 特性

- 支持球面网格和平面网格
- 支持各种地学模式网格
- 支持动态扩展和并行一致性检查
- 支持MPI和OpenMP并行
- 高效并行实现、低计算冗余

## 依赖库及软件

- GNU C++编译器或Intel C++编译器
- MPI
- GUN Make
- Google Test（可选）
- OpenCV（可选）
- NetCDF（可选）

## 编译

1. 根据本机实际情况修改Makefile中的`CXX`变量及`MPI_PATH`变量
2. 在软件目录中执行 `make` 

## 执行

可以使用如下命令运行PatCC：

 `OMP_NUM_THREADS=nt mpiexec -n np ./patcc gridFile`

并相应替换以下参数：

- **nt**：OpenMP线程数
- **np**：MPI进程数
- **GridFile**：符合规定格式的网格

三角化结束后，程序会将三角化结果写入 `log/global_triangles_*` 文件中

## 网格文件格式

PatCC1采用ASCII文本格式作为网格输入文件，具体格式要求如下：

- 第一行：网格总点数
- 第二行：网格边界值，依次为西、东、南、北边界，使用空格分隔
- 第三行及之后的每行代表一个网格点的经纬坐标，使用空格分隔

文件中的坐标值均采用角度为单位，该目录下的`test.grid`为网格示例文件。