作者：王骁

email：wangxiaodm@qq.com

# Monte Carlo
***Q1：计算积分$I=\int^{b}_{a}f(x)dx$。***
计算策略：利用**大数定律**，对某种分布$p(x):\int^{b}_{a}p(x)dx=1$，积分变为： $I=\int^{b}_{a}\frac{f(x)}{p(x)}\cdot p(x)dx=E\left( \frac{f(X)}{p(X)} \right)\approx \frac{1}{n}\sum_{i=1}^{n} \frac{f(x_i)}{p(x_i)}$

***

现在问题转变为：
***Q2：如何获得分布为p(x)的采样？*** 对于简单函数或常见分布，采用反函数法，构造函数法等等； 对于较为复杂的函数，采用接受-拒绝采样法。

但是对于多维分布，或多维函数形式未知时，上面的方法就不太能处理。

***

***Q3：如何对多维分布进行采样？*** 办法就是利用markov chain的收敛性。

# Markov Chain

状态转移矩阵$P_{ij}$表示从状态 i 转移到状态 j 的概率，状态可能是离散的（$i,j=1,2,3,...,n$），也可能时连续的。

若初始处于各个状态的概率为$A_t=(a_1,a_2,...,a_n),\sum\limits_{k=1}^na_k=1$，则下一时刻所处状态的概率分布为：$A_{t+1}=A_tP$。

**Assumption**：状态转移的概率只依赖于它的前一个状态。$P(X_{t+1}|X_t,X_{t-1},...X_0)=P(X_{t+1}|X_t)$。

**Property**：

1.  $\sum_{j=1}^{n}P_{ij}=1$ ；

2.  对于一个非周期的Markov Chain若存在状态转移矩阵$P$且任意两个状态之间连通(**任何两个状态是连通**指的是从任意一个状态可以通过有限步到达其他的任意一个状态，不会出现条件概率一直为0导致不可达的情况。)，则该链可以收敛。$\lim_{n\to\infty}P_{ij}^n=\pi_j$，与初始状态无关；

3.  收敛之后，满足平衡条件：$\pi_j=\sum\limits_{i=1}^{n}{\pi_iP_{ij}},\quad \pi(x_{j})=\int{\pi(x_{i})P(x_{j}|x_{i})dx_{i}}\Rightarrow\pi=\pi P$；

因为Markov Chain最后会收敛到分布$\pi(x)$，所以我们可以取$\pi(x)=p(x)$，只要找到其对应的**状态转移矩阵$P$**，我们就可以从任意初态通过P采样得到我们需要的分布。

***

**Q3.1：如何寻找P？** 若满足细致平衡条件：$\pi(i)P(i,j)=\pi(j)P(j,i)$，则$\pi(x)$是P的平稳分布。

但是一般情况下，细致平衡条件很难达到：$\pi(i)Q(i,j)\neq\pi(j)Q(j,i)$。 令$\alpha(i,j)=\pi(j)Q(j,i)$，则有：

$$
\pi(i)Q(i,j)\alpha(i,j)=\pi(j)Q(j,i)\alpha(j,i)
$$

则分布$\pi(x)$对应的状态转移矩阵为：$P(i,j)\equiv Q(i,j)\alpha(i,j)$。

> 也就是说，我们的目标矩阵 P 可以通过任意一个马尔科夫链状态转移矩阵Q乘以$α(i,j)$得到。$α(i,j)$我们有一般称之为接受率，取值在0\~1之间，可以理解为一个概率值。即目标矩阵 P 可以通过任意一个马尔科夫链状态转移矩阵 Q 以一定的接受率获得。
>
> 这个很像接受-拒绝采样，那里是以一个常用分布通过一定的接受-拒绝概率得到一个非常见分布，这里是以一个常见的马尔科夫链状态转移矩阵 Q 通过一定的接受-拒绝概率得到目标转移矩阵 P ，两者的解决问题思路是类似的。

大致采样步骤为：

1.  首先给定一个任意的状态转移矩阵 Q ，目标分布为$\pi(x)$；

2.  对于任意的一个初态$x_0$，先利用Q来进行一次采样得到下一个状态；

3.  对于第 t 步状态$x_t$，我们利用条件概率分布$Q(x|X=x_t)$得到临时状态$x_*$；

4.  再从均匀分布$U(0,1)$中采样一个值u；

5.  比较u与$\frac{P(x_t,x_*)}{Q(x_t,x_*)}=\alpha(x_t,x_*)=\pi(x_*)Q(x_*,x_t)$之间的大小，

    *   若$u\textless \alpha(x_t,x_*)$，接受，$x_{t+1}=x_*$；
    *   否则，拒绝，$x_{t+1}=t_t$，继续进行第二步。

6.  假设达到平衡状态的采样阈值为$n_1$，总采样次数为$n_2$，则$\{x_{n_1},x_{n_1+1},...x_{n_2}\}$为目标分布的采样。ps：好像没有规定$Q(i,j)\textgreater P(i,j)$？

***

接受率$\alpha(i,j)$可能过小，导致收敛过慢，效率太低。

**Q3.2 如何提高采样效率？**

# sampling method

## M-H采样（Metropolis-Hastings）

在细致平衡条件两侧乘上同一个因子，使得$\max\{C\alpha(i,j),C\alpha(j,i)\}=1$，所以新的接受率为：

$$
\alpha(i,j)=\min\left\{ \frac{\pi(j)Q(j,i)}{\pi(i)Q(i,j)},1 \right\}
$$

一般我们会取Q为对称矩阵：$Q^T=Q$，接受率简化为$\alpha(i,j)=\min\left\{ \frac{\pi(j)}{\pi(i)},1 \right\}$ 如，$Q(i,j)=\frac{1}{\sqrt{2\pi}}e^{-\frac{(x_j-x_i)^2}{2}}$，前后状态$x_i,x_j$也可以是高维的。

***

**Q3.3高维时，采样效率仍然很低，且高维分布$\pi(\vec{x})$常常未知，只知道某个维度的条件概率$P(x_i|x_1,...x_{i-1},x_{(i+1)},...x_n)$。**

## Gibbs采样

Gibbs采样是M-H采样的一个特例。M-H中我们取Q为对称矩阵，这里我们取另一种特殊的Q使得接受率$\alpha=1$。

对于高维分布：$\pi(x_1,x_2,...,x_n)$，假定t时刻处于某个态：$\vec{X}^{(t)}=(x_1^{(t)},x_2^{(t)},...,x_i^{(t)},...x_n^{(t)})$， 我们取$Q(x_i^{(t)},x_i^{(t+1)})=P(x_{i}^{(t)} \to x_i^{(t+1)}|x_1^{(t)},x_2^{(t)},...,x_{i-1}^{(t)},x_{i+1}^{(t)},...,x_n^{(t)})=P(x_{i}=x_i^{(t+1)}|x_{-i}^{(t)})$。 **这里似乎隐含了一个条件，即沿某个轴的转移概率不受上一个状态的影响？** 则有：

$$
\begin{align*}
\pi(x_i^{(t)},x_{-i}^{(t)})Q(x_i^{(t)},x_i^{(t+1)})=P(x_{-i}^{(t)})P(x_i^{(t)}|x_{-i}^{(t)})P(x_i^{(t+1)}|x_{-i}^{(t)})\\
\pi(x_i^{(t+1)},x_{-i}^{(t)})Q(x_i^{(t+1)},x_i^{(t)})=P(x_{-i}^{(t)})P(x_i^{(t+1)}|x_{-i}^{(t)})P(x_i^{(t)}|x_{-i}^{(t)})
\end{align*}
$$

可以看到，上面两式左侧相等，意味着当只在某个轴上改变状态时，一定满足细致平衡条件。 此时接受率$\alpha(i,j)\equiv1$。

**summary:** Gibbs采样的状态转移矩阵可以表示为：

$$
Q(\vec{x}^0,\vec{x}^*)=\left\{
\begin{align*}
P(x_{i}=x_{i}^{*}|x_{-i}=x^{0}_{-i})&,\quad 仅第i个分量变化，其余分量不变:x^{0}_{-i}=x^{*}_{-i}\\
0&,\quad 多余一个分量发生变化
\end{align*}
\right.
$$

    注意：只有将每个轴都采样一次之后得到的才是我们需要的样本，只采样一个轴的Markov Chain是可约的。

***

MH与Gibbs可以结合 ![](Figs/MH+Gibbs.png)


**Q4随着参数维度越来越高，MH的采样方法效率太低，random walk的proposal不再合适，我们需要更加有指向性的采样方法**
# HMC(Harmitonian/Hybrid MC)

文献：

*   [A Conceptual Introduction to Hamiltonian Monte Carlo](zotero://select/library/items/6ISMDDBQ)

*   [MCMC using Hamiltonian dynamics](zotero://select/library/items/X7D39R6B)python库：[pyhmc](https://pythonhosted.org/pyhmc/ "https://pythonhosted.org/pyhmc/")改进M-H在高维时采样效率过慢的问题，我们使用HMC方法代替random walk。

## HMC基本原理

首先明确，$\vec{q}$ 是我们要采样的n维参数矢量，我们的目的是得到参数的某个特定分布$\pi(\vec{q})$ 。

### Hamitonian mechanism

哈密顿力学中，我们有哈密顿方程：

$$
\left\{
\begin{align*}
\frac{\partial{H}}{\partial{p_i}}&=\dot{q}_i\\
\frac{\partial{H}}{\partial{q_i}}&=-\dot{p}_i
\end{align*}
\right.
$$

其中，Hamitonian $H=H(q,p)=K(q)+U(p)$，是位置和动量的函数。

### 扩展到概率空间

对概率引入类似的哈密顿量。

为了利用哈密顿量，我们将参数空间从q扩展到相空间$(q,p)$，相空间的联合概率分布称为**canonical distribution**，$\pi(q,p)=\pi(p|q)\pi(q)$，它不依赖于参数的具体选取形式，只与系统的能量有关，故可以写成：

$$
\begin{align*}
\pi(q,p)&=\frac{1}{Z}e^{-H(q,p)}\\
\Rightarrow H(q,p)&=-\ln{\pi(q,p)}+const\\
&=-\ln\pi(p|q)-\ln\pi(q)+const\\
&\equiv K(p,q)+U(q)+const
\end{align*}
$$

将相空间每个状态记作$(q,p)$，势能$U(q)=-\ln{\pi(q)}$，动能$K(q,p)$的形式有多种，其中最常用的形式是**Euclidean-Gaussian Kinetic Energies** ：

$$
K(p,q)=K(p)=\frac{1}{2}p^TM^{-1}p+\log{|M|}+const
$$

M是**对称、正定**的，而且往往是对角矩阵，物理意义上是质量矩阵，常取单位矩阵乘上标量因子，故可放入常数项中$K(p)=\frac{1}{2}p^TM^{-1}p$。

忽略常数项后，很容易发现对应的条件分布$\pi(p|q)$是*高斯分布*，M对应协方差矩阵，中心在0点。

势能的形式与目标分布$\pi(q)$有关，一般我们取后验分布

$$
\pi(q)=P(q|Data)\propto P(Data|q)P(q)
$$
所以势能：

$$
U(q)=-\ln[P(q)\cdot P(Data|q)]
$$

最后，哈密顿方程可以写成：

$$
\begin{align*}
\dot{\vec{q}}&=M^{-1}\vec{p}\\
\dot{\vec{p}}&=-\frac{\partial{U}}{\partial{\vec{q}}}=-\nabla_\vec{q}U
\end{align*}
$$

***

给予能量守恒条件，$H(q,p)=E$，将正则分布分解为**微正则分布(microcanonical distribution)** 和**边缘能量分布(marginal energy distribution)** 的乘积：$\pi(q,p)=\pi(\theta_E|E)\pi(E)$。

通过固定位置对动量的采样$\pi(p|q)$，我们可以得到能量的条件分布：$\pi(E|q)$，称为**energy transition distribution**。能量转移分布与边缘能量分布的形状应该尽量接近。

### 微分方程求解

微分方程的求解比较困难，常用数值方法求解，如：Euler’s method，A modification of Euler’s method。在HMC中我们使用The leapfrog method。 
	这里是一个可以优化的地方。

t时刻位于：$(q(t),p(t))$，则经过一小段时间$\epsilon$后，到达相空间$(q(t+\epsilon),p(t+\epsilon))$，“蛙跳法”的算法如下：

$$
\begin{align*}
p_i\left( t+\frac{\epsilon}{2} \right)&=p_i(t)-\frac{\epsilon}{2} \frac{\partial{U}}{\partial q_i}(q(t))\\
q_i(t+\epsilon)&=q_i(t)+\epsilon \frac{\partial{K}}{\partial p_i}\left( p\left( t+\frac{\epsilon}{2} \right) \right)\\
p_i(t+\epsilon)&=p_i\left( t+\frac{\epsilon}{2} \right)-\frac{\epsilon}{2} \frac{\partial{U}}{\partial q_i}(q(t+\epsilon))
\end{align*}
$$

这是每次迭代的计算单元，注意到，赏赐迭代的最后一步与下次迭代的第一步能够合并：

$$
\begin{align*}
p_i\left( t+\epsilon+\frac{\epsilon}{2} \right)&=p_i\left( t+\epsilon \right)-\frac{\epsilon}{2} \frac{\partial{U}}{\partial q_i}(q(t+\epsilon))\\
&=p_i\left( t+\frac{\epsilon}{2} \right)-\epsilon \frac{\partial{U}}{\partial q_i}(q(t+\epsilon))
\end{align*}
$$

因此可以将算法简化为除了第一次和最后一次只将$p_i$迭代半步外，中间依次将$q_i,p_i$迭代一个$\epsilon$ 。

若$\pi(p|q)$取高斯形式，则$K(p)=\frac{1}{2}p^TM^{-1}p=\frac{1}{2}p_i(M^{-1})_{ij}p_j,\quad \frac{\partial{K}}{\partial p_k}=(M^{-1})_{kj}p_j$ ，M为动量各分量之间的协方差矩阵。

## 算法及实现

MCMC中有关HMC的部分可以以以下方式引入：

1.  首先从初始点q开始，我们引入辅助动量$p\sim N(\vec{0},M)$，得到相空间的一组初始点：$(q,p)$；
2.  从此点开始，利用HMC，经过L步步长为$\epsilon$的迭代获得 proposal$(q^*,p^*)$，有关proposal的获得按照如下步骤：
	1. 已知：目标分布为后验分布$\pi(\vec{\theta}|Data)$，正则分布为$\pi(\vec{\theta},\vec{p}_\theta)=\pi(\vec{p}_\theta|\vec{\theta})\cdot\pi(\vec{\theta}|Data)$，系统的哈密顿量$H(\vec{\theta},\vec{p}_\theta)=-\ln\pi(\vec{\theta},\vec{p}_\theta)+const$，所以动能项$K(\vec{p}_\theta)=\frac{1}{2}p^TM^{-1}p=\frac{1}{2}p_i(M^{-1})_{ij}p_j$，势能项$U(\vec{\theta})=-\ln[P(\vec{\theta})\cdot P(Data|\vec{\theta})]$ 。
	2. 哈密顿方程为：$$\begin{align*}
\dot{\theta}_i&=M^{-1}_{ij}p_j\\
\dot{p_i}&=-\frac{\partial{U}}{\partial{\theta_i}}
\end{align*}$$利用数值法求解该方程得到下一个点$(\vec{\theta}^*,\vec{p}_\theta^*)$。
3. 由MH采样中的方法，我们以接受率为$\alpha=\min\left\{ 1, \frac{\pi(q^*,p^*)}{\pi{(q,p)}} \right\}=\min\{1, \exp{(-H(q^*,p^*)+H(q,p))}\}$的概率来接受该采样。接受率理论上为1，但是由于数值方法存在误差，可能不为一，但接受率依然很大。



代码实现：[HMC realize](D:/My_All_Code/Code_python/python_learning/emcee/HMC_Move.py)

注：$$\begin{align*}
\alpha&=\exp{(-H(q^*,p^*)+H(q,p))}\\
&=\frac{\pi(q^*)Q^*}{\pi(q)Q}\\
\Rightarrow \frac{Q^*}{Q}&=\alpha \frac{\pi(q)}{\pi(q^*)}\\
\Rightarrow \ln{\frac{Q^*}{Q}}&=\ln{\alpha}+\ln{\pi(q)}-\ln{\pi(q^*)}\\
&=-T^*-U^*+T+U-U+U^*\\
&=T-T^*
\end{align*}$$

# Bayes
一般来说，为了将模型拟合到N个数据上，我们认为似然函数满足以模型理论预测为中心，每个点的误差棒组成的协方差的多维高斯分布：
$$
Likelihood((x_{exp},y_{exp}),\theta)=\frac{1}{\sqrt{2\pi}|\Sigma|}\exp{(-\frac{1}{2}[y_{exp}-f(x_{exp};\theta)]^T\cdot\Sigma^{-1}\cdot[y_{exp}-f(x_{exp};\theta)])}
$$
其中，协方差矩阵可以简单地认为是由数据点的误差棒构成：
$$
\Sigma=\begin{pmatrix}\sigma_1&0&...&0 \\
0&\sigma_2&...&0 \\
...&...&...&... \\
0&0&...&\sigma_n\end{pmatrix}
$$
则似然函数可以简化为：
$$
Likelihood((x_{exp},y_{exp}),\theta)=\frac{1}{\sqrt{2\pi}\prod\limits_i\sigma_i}\exp{\left( -\frac{1}{2}\sum\limits_i\left[ \frac{y_i-f(x_i;\theta)}{\sigma_i} \right]^2 \right)}
$$

