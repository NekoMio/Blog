---
title: "「LibreOJ Round #6」花火"
date: 2018-4-11 10:46:10
categories: 题解
tags: 
    - DP
---

#### 题目描述

>「Hanabi, hanabi……」  
>一听说祭典上没有烟火，Karen 一脸沮丧。  
>「有的哦…… 虽然比不上大型烟花就是了。」  
>还好 Shinobu 早有准备，Alice、Ayaya、Karen、Shinobu、Yoko 五人又能继续愉快地玩耍啦！  
>「噢……！不是有放上天的烟花嘛！」Karen 兴奋地喊道。  
>「啊等等……」Yoko 惊呼。Karen 手持点燃引信的烟花，「嗯？？」  
>Yoko 最希望见到的是排列优美的烟火，当然不会放过这个机会…… 不过时间似乎已经不多了。  

$n$ 个烟火排成一排，从左到右高度分别为 $h\_1, h\_2, \cdots, h\_n$，这些高度两两不同。  
每次 Yoko 可以选择两个相邻的烟火交换，这样的交换可以进行任意多次。  
每次 Yoko 还可以选择两个不相邻的烟火交换，但这样的交换至多进行一次。  
你的任务是帮助 Yoko 用最少次数的交换，使这些烟火从左到右的高度递增。  

#### 输入格式

第一行包含一个正整数 $n$。  
第二行包含 $n$ 个正整数 $h\_1, h\_2, \cdots, h\_n$，相邻整数之间用一个空格隔开。  

#### 输出格式
输出一个整数，表示最少的交换次数。  

#### 样例
  ##### 样例输入
  ```plain
  5
  3 5 4 1 2
  ```
  ##### 样例输出
  ```plain
  5
  ```
  ##### 样例解释
  一开始，$5$ 个烟火的高度依次为 $3,5,4,1,2$。  
  第 $1$ 次，交换第 $4$ 根烟火和第 $5$ 根烟火，交换后烟火的高度依次为 $3,5,4,2,1$。  
  第 $2$ 次，交换第 $3$ 根烟火和第 $4$ 根烟火，交换后烟火的高度依次为 $3,5,2,4,1$。  
  第 $3$ 次，交换第 $1$ 根烟火和第 $2$ 根烟火，交换后烟火的高度依次为 $5,3,2,4,1$。  
  第 $4$ 次，交换第 $2$ 根烟火和第 $3$ 根烟火，交换后烟火的高度依次为 $5,2,3,4,1$。  
  第 $5$ 次，交换第 $1$ 根烟火和第 $5$ 根烟火，交换后烟火的高度依次为 $1,2,3,4,5$。  
  可以证明这是交换次数最少的方案。

#### 数据范围与提示
对于所有数据，满足 $1 \leq n \leq 300,000, 1 \leq h\_i \leq n, h\_i$互不相同。  

|Subtask #|分值|$n$|
|:-:|:-:|:-:|
|1|$6$|$\leq 4$|
|2|$11$|$\leq 8$|
|3|$16$|$\leq 100$|
|4|$8$|$\leq 300$|
|5|$13$|$\leq 700$|
|6|$7$|$\leq 2,000$|
|7|$6$|$\leq 6,000$|
|8|$14$|$\leq 60,000$|
|9|$19$|$\leq 300,000$|

## 题解

### Subtask 1 & 2
直接去搜吧，记忆化一下  
因为与$DP$没关系， 就不写了  

****
### Subtask 3
需要一个多项式复杂度的做法  
我们可以证明  
1. **总是可以假设交换不相邻元素在第一次进行。**  
**证明:** 假设交换的两个值为 $x$ 和 $y$ 那么， 我们可以在第一次就交换他们然后可以构造出一个与原先交换次数相同的解  
2. **如果已经交换了不相邻元素，之后的最小交换次数等于此时排列的逆序对数。**  
**证明:** 设这一次交换的数为 $h[i]$ 与 $h[i + 1]$ 则这一次交换后他们与其他数$(j \neq i, j \neq i + 1)$的逆序对数是不变的。 那么改变的只有他们之间的关系。 又因为一个含有逆序对的数列中一定存在 $h[i] > h[i + 1]$ 所以只有不断的找到这样的数对交换即可， 每次交换减少一个逆序对， 交换次数为逆序对数。

可以 $n^2$ 枚举第一次交换的两个数， 然后 $n^2$ 求逆序对  
时间复杂度$O(n^4)$  
[33pts](https://loj.ac/submission/87860)

****
### Subtask 4
与上面的相同， 只是用随便的一些方法将逆序对用 $n \log(n)$ 求  
[41pts](https://loj.ac/submission/87888)

****
### Subtask 5
这时我们需要一个 $O(n^3)$ 的做法  
因为每次交换所改变的逆序对数很少， 只要先求出来再 $O(n)$ 维护即可  
[54pts](https://loj.ac/submission/87943) 

**** 
### Subtask 6 & 7
每次改变的逆序对并不用 $O(n)$ 求  
将 $(x, h\_x)$ 看做平面上的点  
在第一次交换的时候， 只有当 $h[i] > h[j]$ 时才有意义  
然后我们发现这时减少的逆序对数就是以 $(i, h[i])$ 为左上角 $(j, h[j])$ 为右下角的矩形内部点数的二倍  
然后用一个二维前缀和就可以了  
[67pts](https://loj.ac/submission/88043)

****
### Subtask 8 & 9
到这里其实就比较显然了  
我们想要找的是一对点作为左上角与右下角使得其矩形内的点数最多。  
这显然和我们做过的一道题一样， 然后用线段树维护就好了  
时间复杂度 $O(n\log(n))$  
不知道为什么我的程序又长又慢 2333

```c++
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <set>
using namespace std;
inline int read()
{
    int x=0,f=1;char ch=getchar();
    while (ch<'0'||ch>'9'){if(ch=='-')f=-1;ch=getchar();}
    while (ch>='0'&&ch<='9'){x=x*10+ch-'0';ch=getchar();}
    return x*f;
}
const int MAXN = 300005;
#define size(_) ((_) ? (_)->s : 0)
#define _Max(_) ((_) ? (_)->Max : 0)
struct Node
{
    Node *ch[2];
    int s, v, key;
    int val, Max, lazy;
    Node(int x)
    {
        s = 0, v = x;
        key = rand();
        val = lazy = Max = 1;
        ch[0] = ch[1] = NULL;
    }
    void Pushup()
    {
        s = size(ch[0]) + size(ch[1]) + 1;
        Max = max(max(_Max(ch[0]), _Max(ch[1])), val);
    }
    void Pushlazy(int x)
    {
        val += x;
        lazy += x;
        Max += x;
    }
    void Pushdown()
    {
        if (lazy)
        {
            if (ch[0]) ch[0]->Pushlazy(lazy);
            if (ch[1]) ch[1]->Pushlazy(lazy);
            lazy = 0;
        }
    }
}*root;
Node *Merge(Node *A, Node *B)
{
    if (!A) return B;
    if (!B) return A;
    if (A->key < B->key)
    {
        A->Pushdown();
        A->ch[1] = Merge(A->ch[1], B);
        A->Pushup();
        return A;
    }
    else
    {
        B->Pushdown();
        B->ch[0] = Merge(A, B->ch[0]);
        B->Pushup();
        return B;
    }
}
typedef pair<Node*, Node*> DNode;
DNode Split(Node *rt, int k)
{
    if (!rt) return DNode(NULL, NULL);
    DNode o;
    rt->Pushdown();
    if (size(rt->ch[0]) >= k)
    {
        o = Split(rt->ch[0], k);
        rt->ch[0] = o.second;
        rt->Pushup();
        o.second = rt;
    }
    else
    {
        o = Split(rt->ch[1], k - size(rt->ch[0]) - 1);
        rt->ch[1] = o.first;
        rt->Pushup();
        o.first = rt;
    }
    return o;
}
int Rank(Node *rt, int x)
{
    if (!rt) return 0;
    return x <= rt->v ? Rank(rt->ch[0], x) : Rank(rt->ch[1], x) + size(rt->ch[0]) + 1;
}
void Insert(int x)
{
    int k = Rank(root, x);
    DNode y = Split(root, k);
    Node *n = new Node(x);
    root = Merge(y.first, Merge(n, y.second));
}
void Update(int x, int c)
{
    int k = Rank(root, x);
    DNode y = Split(root, k);
    if (y.second) y.second->Pushlazy(c);
    root = Merge(y.first, y.second);
}
void Update(int l, int r, int c)
{
    int k = Rank(root, l);
    DNode x = Split(root, k);
    int k2 = Rank(x.second, r + 1);
    DNode y = Split(x.second, k2);
    if (y.first) y.first->Pushlazy(c);
    root = Merge(x.first, Merge(y.first, y.second));
}
int Calc_Max(Node *rt)
{
    if (!rt->ch[1]) return rt->v;
    else return Calc_Max(rt->ch[1]);
}
int Query(int x)
{
    int k = Rank(root, x);
    DNode y = Split(root, k);
    int ans = _Max(y.second);
    root = Merge(y.first, y.second);
    return ans + 1;
}
set<int> t1, t2;
struct data
{
    int v, id, Mx;
    bool operator < (const data &a) const
    {
        return v < a.v;
    }
};
set<data> t3;
int n;
int a[MAXN], b[MAXN];
long long S[MAXN];
#define lowbit(_) ((_) & (-_))
void add(int x, int c)
{
    for (int i = x; i <= n; i += lowbit(i))
        S[i] += c;
}
long long Q(int x)
{
    long long ans = 0;
    for (int i = x; i >= 1; i -= lowbit(i))
        ans += S[i];
    return ans;
}
int st1[MAXN], st2[MAXN], top1, top2;
int main()
{
    // freopen ("columns3_2.in", "r", stdin);
    // freopen ("columns.out", "w", stdout);
    n = read();
    for (int i = 1; i <= n; i++)
        b[i] = a[i] = read();
    sort(b + 1, b + n + 1);
    for (int i = 1; i <= n; i++)
        a[i] = lower_bound(b + 1, b + n + 1, a[i]) - b;
    long long Sum = 0;
    for (int i = 1; i <= n; i++)
        Sum += Q(n) - Q(a[i]), add(a[i], 1);
    st1[++top1] = 1;
    t1.insert(1);
    for (int i = 2; i <= n; i++)
        if (a[i] > a[st1[top1]])
            st1[++top1] = i, t1.insert(i);
    st2[++top2] = n;
    t2.insert(n);
    for (int i = n - 1; i >= 1; i--)
        if (a[i] < a[st2[top2]])
            st2[++top2] = i, t2.insert(i);
    reverse(st2 + 1, st2 + top2 + 1);
    int ans = 0;
    for (int i = 1; i <= n; i++)
    {
        if (!t1.count(i) && !t2.count(i))
        {
            Update(a[i], 1);
            t3.insert((data){a[i], i, Calc_Max(root)});
        }
        else
        {
            if (t1.count(i))
                Insert(a[i]);
            if (t2.count(i))
            {
                for (auto x : t3)
                    if (x.v > a[i]) break;
                    else
                        Update(x.v, x.Mx, -1);
                auto x2 = t3.begin();
                if (x2 != t3.end()) x2++;
                for (auto x = t3.begin(); x != t3.end(); )
                {
                    if (x->v > a[i]) break;
                    else
                        t3.erase(x);
                    x = x2;
                    if (x2 != t3.end())
                        x2++;
                }
                ans = max(ans, Query(a[i]));
            }
        }
    }
    printf ("%lld\n", min(Sum, Sum - (ans - 2) * 2));
}
```