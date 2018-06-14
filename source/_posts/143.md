---
title: "「LibreOJ Round #8」MINIM"
date: 2018-4-11 20:27:50
categories: 题解
tags: 
    - DP
---

#### 题目描述
取石子游戏的规则是这样的：有若干堆石子，两个玩家轮流操作，每个玩家每次要选一堆取走任意多个石子，但不能不取，无石子可取者输。

现在共有 $n$ 堆石子，其中第 $i$ 堆的数量为 $l\_i$，现在 LCR 需要在每一堆中扔掉一部分（可以不扔也可以全扔），如果第 $i$ 堆的石子在 LCR 操作后**还有剩余**，LCR 就需要付出 $v\_i$ 的代价。LCR 操作完成后神犇会搬来新的一堆个数在 $[0,m]$ 之间的石子，两人玩取石子游戏，LCR 先手。神犇搬运新的一堆石子时会保证自己（后手）必胜，如果他无法做到这一点，就会立即结束游戏。

现在 LCR 有 $q$ 次询问，每次给出一个 $c\in [0,m]$，请你回答如果要让神犇搬来的石子数为 $c$（不能让神犇结束游戏，即使这里要求 $c=0$），LCR 付出代价的总和至少是多少。如果 LCR 不可能通过调整石子使得神犇搬来的石子数为 $c$，输出 `-1`。

NIM is a game of strategy in which two players take turns removing stones from distinct piles. On each turn, a player must remove at least one stone, and may remove any number of stones provided they all come from the same pile. **The player who has no stones to remove loses.**

There are $n$ piles of stones, the $i$-th pile has $l\_i$ stones. Alice needs to remove some of the stones from each pile (removing zero or all of the stones from a certain pile is allowed). If the $i$-th pile remains at least one stone after this operation, Alice has to pay a price of $v\_i$. After Alice's operation, Bob will create a pile of stones whose number is in $[0,m]$ and add it to the game to make sure that Alice — the player who moves first in this NIM game will lose.If he can't ensure that Alice will lose,he will exit from the game immediately.

Now Alice has $q$ queries. Each query gives an integer $c$, and you need to calculate the minimal total price Alice needs to pay to ensure Bob’s new pile has exactly $c$ stones(making Bob exit from the game isn't allowed,even if $c=0$). If this is impossible, print `-1` instead.

#### 输入格式
第一行两个正整数 $n,m$。  
接下来 $n$ 行每行两个正整数 $v\_i,l\_i$ 表示该堆石子的代价和数量。  
接下来一行一个正整数 $q$。  
接下来 $q$ 行每行一个正整数 $c$ 表示询问。  

The first line contains two integers $n,m$.  
Each of the following $n$ lines contains two integers $v\_i,l\_i$.  
The next line contains an integer $q$. Each of the following $n$ lines contains an integers $c$.  

#### 输出格式
输出共 $q$ 行，依次表示每次询问的答案，无解输出 `-1`。

Output contains $q$ lines,containing the answer of each query in order.  

#### 样例
##### 样例输入 Sample Input
```
4 6
2 3
4 4
3 5
5 2
7
0
1
2
3
4
5
6
```

##### 样例输出 Sample Output
```
0
2
2
2
3
3
5
```

#### 数据范围与提示

对于所有数据，$1 \leq n,q \leq 10^5, 1 \leq v\_i \leq 10^9, 0 \leq l\_i \leq m \leq 10^9, c \in [0, m]$。  

详细的数据限制及约定如下（留空表示和上述所有数据的约定相同）

|Subtask #|分值（百分比）|$n,q$|$l\_i,m$|特殊性质|
|:-:|:-:|:-:|:-:|:-:|
|1|15|$\leq 10$|$\leq 10$|-|
|2|20|$\leq 100$|$\leq 100$|-|
|3|20|-|-|对于每个 $i$ 存在非负整数 $k$ 满足 $l\_i=2^k-1$|
|4|20|$\leq 20000$|$\leq 20000$|$l\_i,v\_i$ 在范围内均匀随机（使用 `std::mt19937` 并对最大值取模）|
|5|25|-|-|$l\_i,v\_i$ 在范围内均匀随机（使用 `std::mt19937` 并对最大值取模）|


## 题解
### Subtask 1
直接爆搜  

### Subtask 2

考虑 $DP$  
因为想要使得后手必胜，则初始状态的异或值为 $0$ , 那么就可以 $DP$ 了, 如果想要让神犇搬来的石子数为 $c$ 则前面的异或值需要为 $c$  
然后定义 $f[i][j]$ 表示 $DP$ 到 $i$ 异或值为 $j$ 的最小花费， 暴力转移即可  
时间复杂度$O(nm^2)$  
[35pts](https://loj.ac/submission/88276)

### Subtask 3

这部分所有的 $i$ 存在非负整数 $k$ 满足 $l\_i=2^k-1$   
可以发现最优解一定是选一个  
我们只要找大于等于 $c$ 的数中 $v\_i$ 的最小值就可以了  
[55pts](https://loj.ac/submission/88286) 

### Subtask 4 & 5

我们发现 [Subtask 2](#Subtask-2) 中的每一层 $f$ 数组是分段的  
当数据随机时可以证明: **段数的期望是常数**  
[CommonAnts的证明](https://loj.ac/article/322)是这样的
>考虑最小的一对满足 $m$ 以内二进制最高位均为 $1$ 的 $l\_i,l\_j$，那么这一对可以保证 $DP$ 数组的任意位置不超过 $v\_i+v\_j$。又根据算法四的分析，每一个权值和不超过 $v\_i+v\_j$ 的子集的贡献不超过一段。  
>由于随机，$v\_i$ 的分布是均匀的，故枚举 $i,j$ 计算概率可得期望段数不超过 $\sum\_{i \leq j}{2^{-j}\sum\_{k=1}^{i+j}{S(i+j-\binom{k-1}{2}, k)}} = O(1)$ ($S$ 表示第二类斯特林数)  
>对于 $m$ 不是 $2$ 的幂的情况，可以分最高位和次高位分别分析，每一位产生的贡献都不超过上述的常数。

然后就可以记录段数转移了  
时间复杂度 $O(n+q)$ (Subtask 4 & 5) 或 $O(nlog(n))$ (Subtask 3) 或 $O(nm)$ (Subtask 1 & 2) 

```c++
#include <cstdio>
#include <cstring>
#include <algorithm>
using namespace std;
inline int read()
{
    int x=0,f=1;char ch=getchar();
    while (ch<'0'||ch>'9'){if(ch=='-')f=-1;ch=getchar();}
    while (ch>='0'&&ch<='9'){x=x*10+ch-'0';ch=getchar();}
    return x*f;
}
inline void gmin(int &a, int b) { if (a > b) a = b; }
const int MAXN = 100005;
int M;
struct data
{
    int v, l;
    data(int _v = 0, int _l = 0) : v(_v), l(_l) {}
    data operator + (const data &a) const 
    {
        data c(v + a.v, l | a.l);
        for (int i = M; i >= 0; i--)
            if (l & a.l & (1 << i))
            {
                c.l |= ((1 << (i + 1)) - 1);
                break;
            }
        return c;
    }
    bool operator < (const data &a) const 
    {
        return v < a.v;
    }
}a[MAXN], s[MAXN], t[MAXN], tmp[MAXN];
int cnt, len;
int Calc(int x)
{
    for (int i = 1; i <= cnt; i++)
        if (s[i].l >= x)
            return s[i].v;
    return -1;
}
int main()
{
    // freopen ("1.out", "w", stdout);
    int n = read(), m = read();
    for (int i = 1; i <= n; i++)
        a[i].v = read(), a[i].l = read();
    for (M = 31; M && (1 << (M - 1)) >= m; M--);
    s[++cnt] = data(0, 0);
    tmp[0] = data(0, -1);
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= cnt; j++) 
            t[j] = s[j] + a[i];
        len = 0;
        int j = 0, k = 0, l;
        while (j < cnt || k < cnt)
        {
            data tp = j >= cnt ? t[++k] : (k >= cnt ? s[++j] : s[j + 1] < t[k + 1] ? s[++j] : t[++k]);
            for (l = 1; l <= len; l++)
                if (tmp[l].l >= tp.l)
                    break;
            if (l == len + 1)
                tmp[++len] = tp;
        }
        for (int j = 1; j <= len; j++)
            s[j] = tmp[j];
        cnt = len;
    }
    int q = read();
    while (q--)
    {
        printf ("%d\n", Calc(read()));
    }
}
```