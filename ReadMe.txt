道德效价在知觉匹配任务中的非对称性
Moral valence in perceptual matching

一、当前范围
本项目目前仅整理实验1a。每个实验单独存放方案、材料、数据说明和分析。
以后实验1b、1c、2a、2b分别使用3_Exp1b、4_Exp1c、5_Exp2a、6_Exp2b。
本次不创建这四个实验目录。7_Reports统一存放研究报告和手稿，目前不设置PPT目录。

二、目录结构
moral valence/
|-- 1_General_Protocol/
|   `-- About_Protocol.txt
|-- 2_Exp1a/
|   |-- 2_1_Study_specific_protocol/
|   |   `-- About_Study_Protocol.txt
|   |-- 2_2_MaterialProc/
|   |   |-- About_MaterialProc.txt
|   |   `-- 2_2_1_Procedure/
|   |-- 2_3_RawData/
|   |   `-- About_RawData.txt
|   `-- 2_4_Analyses/
|       |-- About_Analyses.txt
|       |-- 2_4_1_Data/
|       |-- 2_4_2_Script/
|       |   `-- 1a_full_code_1030.R
|       |-- 2_4_3_Results/
|       |-- 2_4_4_Figures/
|       `-- 2_4_5_ModelCache/  [本地模型缓存，不提交至Git]
|-- 7_Reports/
|   |-- About_Reports.txt
|   `-- 7_1_Manuscript/
|-- moral valence.Rproj
|-- ReadMe.txt
`-- .gitignore

三、代码与模型
主代码：2_Exp1a/2_4_Analyses/2_4_2_Script/1a_full_code_1030.R。
RT使用lognormal分布；SDT使用Bernoulli分布与probit链接。
本次只整理目录和路径，未修改模型公式、先验、筛选标准、抽样设置或统计结论。
所有MCMC使用CmdStanR。桥采样的add_rstan_model仅提供概率密度接口，不进行RStan MCMC。
正文后验差值图为匹配RT三组与dprime三组；补充图增加不匹配RT三组，共九组。
后验差值图在dprime绘制后统一保存。CSV结果表仍按分析小节保存。

四、如何在RStudio运行
1. 先打开本项目根目录的moral valence.Rproj。
2. 再打开主代码。curDir <- getwd()以该R项目为根目录。
3. data_dir仍指向原始数据原位置；见2_Exp1a/2_3_RawData/About_RawData.txt。
4. CmdStan和临时文件配置仍使用D:/Stan_BF；其他电脑需自行配置相应路径与R环境。
5. 需要R、CmdStan及readr、dplyr、tidyr、ggplot2、ggdist、patchwork、brms、cmdstanr、posterior、bayesplot、bridgesampling。
6. 模型保存目录为2_Exp1a/2_4_Analyses/2_4_5_ModelCache。本次没有复制或新拟合模型。
   该目录没有对应缓存时，完整Source会运行MCMC；请按需要运行小节。
7. 普通R启动时的临时目录由用户的.Renviron设置；仅修改工作目录不会改变tempdir()。
   本机分析临时目录应为D:/Stan_BF/R_tmp。

五、现有文件状态
本轮放入1a代码与目录说明。未从其他位置复制原始数据、既有模型、结果、图或Word手稿。
方案、材料、原始数据和手稿目录目前是预留位置；ReadMe说明不代表相关研究文件已归档。
ReadMe中不声明未核实的预注册、伦理批准或开放许可证。
.gitignore排除模型缓存、R工作区和临时文件；结果与原始数据是否分享需另行选择。

六、GitHub Desktop
桌面moral valence是本地整理目录，其中moral valence & perceptual matching是现有GitHub仓库。
本轮已将同一结构和同一版本代码同步到该仓库，保留.git与远程地址。
GitHub Desktop中查看Changes，确认后Commit to main，再Push origin。
本轮没有提交或推送。后续上传以仓库目录内的代码为准；两个文件夹不会自动同步。
