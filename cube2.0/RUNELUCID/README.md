# 运行以ELUCID重构初始密度场为IC的模拟：
## 设置参数
### 在[parameters.f90](./parameters.f90)中设置基础参数
- 格点ng=500，盒子大小box=500Mpc/h
- ratio_cs：粗格点/标准格点长度  2 or 5，在PM2时2或5结果差不多，PP时会有差别
- nnt：每个node中tile数，每维度  2 or 5
- nns：每个tile中subtile数，每维度  4 or 5
- WMAP5 cosmology：h0=0.72，omega_cdm=0.214，omega_bar=0.044，s8=0.807332993

### 在[ic.f90](./utilities/ic.f90)中设置
- 替换随机初始条件为ELUCID初始密度场cxyz_251_500_500.bin 
    - ELUCID_small，替换ng=500的场（一比一替换）
    - ELUCID_big，替换ng>500的场（在傅里叶空间替换中心，等效于小尺度的k是随机场）
- gradphi的计算要用旧的两点差分，新的差分会使密度场整体+0.5，和观测坐标对不上

### 在[kick.f90](./kick.f90)中设置
- PM2: 是否计算PM2的力
- PM3: 是否计算PM3的力
- PP: 是否计算PP力
