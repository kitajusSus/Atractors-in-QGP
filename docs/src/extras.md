# Extras
Sometimes I can find techniques that see usefull for my work, I will save them  here


## Dynamical Systems Analysis

Beyond PCA, AttractorsQGP.jl provides tools to analyze the system from the perspective of chaos theory and dynamical systems.

### Potencial Way for upgrades
While I am doing my research sometimes I find out about interesting algorithms:

#### Lyapunov Exponent

Measures the rate of separation of infinitesimally close trajectories. A negative LLE strongly
indicates the presence of an attracting manifold.

```julia
u0 = [2.0, 5.0] # [T0 in fm⁻¹, A0]
lle = run_LLE(model_brsss, u0, (0.22, 5.0); perturbation=1e-6)
println("Local Lyapunov Exponent: ", lle)

```

#### Intrinsic Dimension

Estimates the intrinsic geometric dimension of the data cloud at a given time
using the participation ratio of the covariance matrix eigenvalues.

```julia
_, X_tau = get_tau_slice(dataset, 0.5)
dim = estimate_dimension(X_tau)

```




## youtube


### Dimensionality

- [The Curse of Dimensionality youtube ](https://www.youtube.com/watch?v=9Tf-_mJhOkU)

How distances Increse in higher dimensions?
How to indentify which features are more significant/important/relevant


### LLE
- [Localy Linear Embeding Lecture](https://www.youtube.com/watch?v=scMntW3s-Wk&t=59s)

Shows how PCA and PCA-kernel got a little problems with holding the structure of data

> How does it work?
> [!IMPORTANT]
> $\epsilon(W) = \sum_i = |\vec{X}_i - \sum_j W_{ij}\vec{Y}_j|^2$
Computing a set of weights that can be used for reconstructing a point










# 

```julia

julia> println("  τ [fm/c]  | Liczba punktów | Efektywny Wymiar (Participation Ratio)")
  τ [fm/c]  | Liczba punktów | Efektywny Wymiar (Participation Ratio)

julia> println("="^65)
=================================================================

julia> for tau in wysk_tau
           _, X_tau = get_tau_slice(dane, tau; feature_cols=[2, 3])

           if size(X_tau, 1) > 2
               X_norm = normalizuj_max(X_tau)

               d_eff = estimate_dimension(X_norm)

               @printf("    %.3f    |      %5d     |                %.4f\n", tau, size(X_tau, 1), d_eff)
           end
       end
    0.200    |       5000     |                1.8816
    0.210    |       5000     |                1.9355
    0.220    |       5000     |                1.9285
    0.230    |       5000     |                1.8682
    0.240    |       5000     |                1.7774
    0.250    |       5000     |                1.6772
    0.260    |       5000     |                1.5822
    0.270    |       5000     |                1.4987
    0.280    |       5000     |                1.4282
    0.290    |       5000     |                1.3698
    0.300    |       5000     |                1.3220
    0.310    |       5000     |                1.2829
    0.320    |       5000     |                1.2509
    0.330    |       5000     |                1.2248
    0.340    |       5000     |                1.2034
    0.350    |       5000     |                1.1858
    0.360    |       5000     |                1.1714
    0.370    |       5000     |                1.1595
    0.380    |       5000     |                1.1496
    0.390    |       5000     |                1.1416
    0.400    |       5000     |                1.1349
    0.410    |       5000     |                1.1294
    0.420    |       5000     |                1.1249
    0.430    |       5000     |                1.1213
    0.440    |       5000     |                1.1183
    0.450    |       5000     |                1.1159
    0.460    |       5000     |                1.1140
    0.470    |       5000     |                1.1126
    0.480    |       5000     |                1.1114
    0.490    |       5000     |                1.1106
    0.500    |       5000     |                1.1100
    0.510    |       5000     |                1.1096
    0.520    |       5000     |                1.1094
    0.530    |       5000     |                1.1093
    0.540    |       5000     |                1.1094
    0.550    |       5000     |                1.1096
    0.560    |       5000     |                1.1099
    0.570    |       5000     |                1.1102
    0.580    |       5000     |                1.1106
    0.590    |       5000     |                1.1110
    0.600    |       5000     |                1.1115
    0.610    |       5000     |                1.1120
    0.620    |       5000     |                1.1125
    0.630    |       5000     |                1.1131
    0.640    |       5000     |                1.1136
    0.650    |       5000     |                1.1142
    0.660    |       5000     |                1.1148
    0.670    |       5000     |                1.1154
    0.680    |       5000     |                1.1160
    0.690    |       5000     |                1.1166
    0.700    |       5000     |                1.1172
    0.710    |       5000     |                1.1177
    0.720    |       5000     |                1.1183
    0.730    |       5000     |                1.1189
    0.740    |       5000     |                1.1195
    0.750    |       5000     |                1.1200
    0.760    |       5000     |                1.1206
    0.770    |       5000     |                1.1211
    0.780    |       5000     |                1.1217
    0.790    |       5000     |                1.1222
    0.800    |       5000     |                1.1227
    0.810    |       5000     |                1.1232
    0.820    |       5000     |                1.1237
    0.830    |       5000     |                1.1242
    0.840    |       5000     |                1.1246
    0.850    |       5000     |                1.1251
    0.860    |       5000     |                1.1255
    0.870    |       5000     |                1.1258
    0.880    |       5000     |                1.1262
    0.890    |       5000     |                1.1265
    0.900    |       5000     |                1.1268
    0.910    |       5000     |                1.1272
    0.920    |       5000     |                1.1275
    0.930    |       5000     |                1.1278
    0.940    |       5000     |                1.1281
    0.950    |       5000     |                1.1283
    0.960    |       5000     |                1.1285
    0.970    |       5000     |                1.1287
    0.980    |       5000     |                1.1289
    0.990    |       5000     |                1.1291
    1.000    |       5000     |                1.1292
    1.010    |       5000     |                1.1294
    1.020    |       5000     |                1.1295
    1.030    |       5000     |                1.1297
    1.040    |       5000     |                1.1298
    1.050    |       5000     |                1.1299
    1.060    |       5000     |                1.1300
    1.070    |       5000     |                1.1300
    1.080    |       5000     |                1.1301
    1.090    |       5000     |                1.1301
    1.100    |       5000     |                1.1302
    1.110    |       5000     |                1.1302
    1.120    |       5000     |                1.1302
    1.130    |       5000     |                1.1303
    1.140    |       5000     |                1.1303
    1.150    |       5000     |                1.1303
    1.160    |       5000     |                1.1303
    1.170    |       5000     |                1.1303
    1.180    |       5000     |                1.1302
    1.190    |       5000     |                1.1302
    1.200    |       5000     |                1.1302
    1.210    |       5000     |                1.1302
    1.220    |       5000     |                1.1302
    1.230    |       5000     |                1.1301
    1.240    |       5000     |                1.1301
    1.250    |       5000     |                1.1301
    1.260    |       5000     |                1.1300
    1.270    |       5000     |                1.1300
    1.280    |       5000     |                1.1300
    1.290    |       5000     |                1.1299
    1.300    |       5000     |                1.1299
    1.310    |       5000     |                1.1298
    1.320    |       5000     |                1.1298
    1.330    |       5000     |                1.1297
    1.340    |       5000     |                1.1297
    1.350    |       5000     |                1.1297
    1.360    |       5000     |                1.1296
    1.370    |       5000     |                1.1295
    1.380    |       5000     |                1.1295
    1.390    |       5000     |                1.1294
    1.400    |       5000     |                1.1294
    1.410    |       5000     |                1.1293
    1.420    |       5000     |                1.1293
    1.430    |       5000     |                1.1292
    1.440    |       5000     |                1.1292
    1.450    |       5000     |                1.1291
    1.460    |       5000     |                1.1290
    1.470    |       5000     |                1.1290
    1.480    |       5000     |                1.1289
    1.490    |       5000     |                1.1289
    1.500    |       5000     |                1.1288
```
