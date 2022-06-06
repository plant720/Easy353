# easy353 更迭版本记录 

## 2022.04.25 1.0.0
基本功能实现：用于NGS二代测序数据中提取特定的序列(主要是被子植物353) \
    包括从raw/clean data中过滤与目标序列相关的read,并将结果输出到指定的文件夹中 \
    将过滤后的read进行assemble,并将结果输出到指定的文件夹中

## 2022.04.25 1.0.1
fix bug:  修改要求输入read的长度 \
fix bug: 修改要求输入read的长度

## 2022.5.1 1.0.2
function add: 增加scaffold

## 2022.5.6 1.0.3
fix bug: filter过滤中设定了read反向互补，但实际上read并未反向互补

## 2022.5.7 1.0.4
fix bug: assemble报错

## 2022.5.7 1.0.5
fix bug: 修改输出文件结构
