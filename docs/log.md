# easy353 更迭版本记录 
> 版本号： a,软件整体架构变动 b,功能新增 c,bug修复 d,github提交次数
## 2022.04.25 1.0.0
基本功能实现：用于NGS二代测序数据中提取特定的序列(主要是被子植物353) \
    包括从raw/clean data中过滤与目标序列相关的read,并将结果输出到指定的文件夹中 \
    将过滤后的read进行assemble,并将结果输出到指定的文件夹中

## 2022.04.25 1.0.1
fix bug:  修改要求输入read的长度

## 2022.5.1 1.1.0
function add: 增加scaffold

## 2022.5.6 1.1.1
fix bug: filter过滤中设定了read反向互补，但实际上read并未反向互补

## 2022.5.7 1.1.2
fix bug: assemble报错

## 2022.5.7 1.1.3
fix bug: 修改输出文件结构

## 2022.6.7 1.2.0
function add: 增加输出结果的统计信息
fix bug: 解决组装结果报错的问题