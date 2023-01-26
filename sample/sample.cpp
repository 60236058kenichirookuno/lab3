#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <lbfgs.h>
#include <GL/glut.h>
#include <vector>
#include <random>
using namespace std;
#define SEGMENT 10 // 一本の糸のセグメント数
#define YARN 20 //糸の総数
#define N   SEGMENT*YARN*2 //全要素（高さと傾き）
#define DIMENSION 10 //分割数
#define SPLINE SEGMENT*YARN*DIMENSION //スプライン曲線全要素
#define RADIUS 0.5 //糸の半径
#define INTERVAL RADIUS*2 //(仮) //同じ方向の糸同士の間隔
#define Ks 2 //
#define Kc 200 //
#define CIRCLE 10 //回転の分割数
#define CIRCLEF 10 //F回転の分割数
#define RATIO INTERVAL*YARN/(2*SEGMENT)  //比率
#define PLY 3 //より糸の本数
#define Rply 0.3 //より糸と糸の距離
#define Aply 5.0 //より糸の螺旋のピッチ
#define Tply 0 //より糸の螺旋の初期位相
#define YPdistance 0.25 //より糸の半径
#define FIBER 10 //繊維の本数
#define Fradius 0.05 //繊維の半径
#define Afiber 0.02 //繊維の螺旋のピッチ
lbfgsfloatval_t K[N];
static const double pi = 3.141592653589793;
lbfgsfloatval_t radius = (lbfgsfloatval_t)RADIUS;
lbfgsfloatval_t circle = (lbfgsfloatval_t)CIRCLE;
lbfgsfloatval_t circlef = (lbfgsfloatval_t)CIRCLEF;
lbfgsfloatval_t F[FIBER][2];
GLfloat o[10][2];
GLfloat O[10][20][2];

int w = 0;
int z = 0;
class objective_function
{
protected:
    lbfgsfloatval_t *m_x;

public:
    objective_function() : m_x(NULL)
    {
    }

    virtual ~objective_function()
    {
        if (m_x != NULL) {
            lbfgs_free(m_x);
            m_x = NULL;
        }
    }

    int run(int n)
    {
        lbfgsfloatval_t fx;
        lbfgsfloatval_t *m_x = lbfgs_malloc(n);

        if (m_x == NULL) {
            printf("ERROR: Failed to allocate a memory block for variables.\n");
            return 1;
        }

        /* Initialize the variables. */
        lbfgsfloatval_t takasaI[YARN / 2][SEGMENT]=
        
        { { 1.0,1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0,-1.0,1.0 },

        { 1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0,-1.0 },

        { -1.0,1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0 },

        { -1.0,-1.0,1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,1.0 },

        { 1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0,-1.0,1.0,1.0 },

        { 1.0,1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0,-1.0,1.0 },

        { 1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0,-1.0 },

        { -1.0,1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0 },

        { -1.0,-1.0,1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,1.0 },

        { 1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0,-1.0,1.0,1.0 } };

        lbfgsfloatval_t takasaJ[YARN / 2][SEGMENT] = 
        
        { { -1.0,-1.0,1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0 },

        { -1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0 },

        { 1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,-1.0,1.0 },

        { 1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,-1.0 },

        { -1.0,1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0 },

        { -1.0,-1.0,1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0 },

        { -1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0 },

        { 1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,-1.0,1.0 },

        { 1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,-1.0 },

        { -1.0,1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0 } };

        lbfgsfloatval_t katamukiI[YARN / 2][SEGMENT] = 
        
        { { 0.0,-2.0,0.0,2.0,0.0,0.0,-2.0,0.0,2.0,0.0 },

        { 0.0,0.0,-2.0,0.0,2.0,0.0,0.0,-2.0,0.0,2.0 },

        { 2.0,0.0,0.0,-2.0,0.0,2.0,0.0,0.0,-2.0,0.0 },

        { 0.0,2.0,0.0,0.0,-2.0,0.0,2.0,0.0,0.0,-2.0 },

        { -2.0,0.0,2.0,0.0,0.0,-2.0,0.0,2.0,0.0,0.0 },

        { 0.0,-2.0,0.0,2.0,0.0,0.0,-2.0,0.0,2.0,0.0 },

        { 0.0,0.0,-2.0,0.0,2.0,0.0,0.0,-2.0,0.0,2.0 },

        { 2.0,0.0,0.0,-2.0,0.0,2.0,0.0,0.0,-2.0,0.0 },

        { 0.0,2.0,0.0,0.0,-2.0,0.0,2.0,0.0,0.0,-2.0 },

        { -2.0,0.0,2.0,0.0,0.0,-2.0,0.0,2.0,0.0,0.0 } };

        lbfgsfloatval_t katamukiJ[YARN / 2][SEGMENT] = 
        
        { { 0.0,2.0,0.0,-2.0,0.0,0.0,2.0,0.0,-2.0,0.0 },

        { 0.0,0.0,2.0,0.0,-2.0,0.0,0.0,2.0,0.0,-2.0 },

        { -2.0,0.0,0.0,2.0,0.0,-2.0,0.0,0.0,2.0,0.0 },

        { 0.0,-2.0,0.0,0.0,2.0,0.0,-2.0,0.0,0.0,2.0 },

        { 2.0,0.0,-2.0,0.0,0.0,2.0,0.0,-2.0,0.0,0.0 },

        { 0.0,2.0,0.0,-2.0,0.0,0.0,2.0,0.0,-2.0,0.0 },

        { 0.0,0.0,2.0,0.0,-2.0,0.0,0.0,2.0,0.0,-2.0 },

        { -2.0,0.0,0.0,2.0,0.0,-2.0,0.0,0.0,2.0,0.0 },

        { 0.0,-2.0,0.0,0.0,2.0,0.0,-2.0,0.0,0.0,2.0 },

        { 2.0,0.0,-2.0,0.0,0.0,2.0,0.0,-2.0,0.0,0.0 } };

        
        int l = 0;
        int j = 0;
        for (int k = 0; k < YARN / 2; k++) {
            j = 0;
            for (int i = l; i < l + 2 * SEGMENT; i += 2) {
                m_x[i] = takasaI[k][j];
                m_x[i + 1] = katamukiI[k][j];

                j = j + 1;
            }
            l += 2 * SEGMENT;
            j = 0;
            for (int i = l; i < l + 2 * SEGMENT; i += 2) {
                m_x[i] = takasaJ[k][j];
                m_x[i + 1] = katamukiJ[k][j];

                j = j + 1;
            }
            l += 2 * SEGMENT;
        }
        /*
            Start the L-BFGS optimization; this will invoke the callback functions
            evaluate() and progress() when necessary.
         */
        int ret = lbfgs(n, m_x, &fx, _evaluate, _progress, this, NULL);

        /* Report the result. */
        //printf("L-BFGS optimization terminated with status code = %d\n", ret);
        //printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, m_x[0], m_x[1]);

        for (int i = 0; i < n; i = i + 1) {
            K[i] = m_x[i];
        }
        
        return ret;
    }

protected:
    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        return reinterpret_cast<objective_function*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        lbfgsfloatval_t fx = 0.0;
        lbfgsfloatval_t d = 0.0;
        lbfgsfloatval_t t0 = 0.0;
        lbfgsfloatval_t t1 = 0.0;
        lbfgsfloatval_t t2 = 0.0;
        lbfgsfloatval_t t3 = 0.0;
        lbfgsfloatval_t t4 = 0.0;
        
        int syoutotu[YARN * YARN / 4][2] =
        { { 0,20 }, { 2,60 }, { 4,100 } ,{ 6,140 }, { 8,180 }, { 10,220 }, { 12,260 }, { 14,300 }, { 16,340 }, { 18,380 } ,
          { 40,22 }, { 42,62 }, { 44,102 }, { 46,142 }, { 48,182 }, { 50,222 }, { 52,262 }, { 54,302 }, { 56,342 }, { 58,382 },
          { 80,24 }, { 82,64 }, { 84,104 }, { 86,144 }, { 88,184 }, { 90,224 }, { 92,264 }, { 94,304 }, { 96,344 }, { 98,384 },
          { 120,26 }, { 122,66 }, { 124,106 }, { 126,146 }, { 128,186 }, { 130,226 }, { 132,266 }, { 134,306 }, { 136,346 }, { 138,386 },
          { 160,28 }, { 162,68 }, { 164,108 }, { 166,148 }, { 168,188 }, { 170,228 }, { 172,268 }, { 174,308 }, { 176,348 }, { 178,388 },
          { 200,30 }, { 202,70 }, { 204,110 }, { 206,150 }, { 208,190 }, { 210,230 }, { 212,270 }, { 214,310 }, { 216,350 }, { 218,390 },
          { 240,32 }, { 242,72 }, { 244,112 }, { 246,152 }, { 248,192 }, { 250,232 }, { 252,272 }, { 254,312 }, { 256,352 }, { 258,392 },
          { 280,34 }, { 282,74 }, { 284,114 }, { 286,154 }, { 288,194 }, { 290,234 }, { 292,274 }, { 294,314 }, { 296,354 }, { 298,394 },
          { 320,36 }, { 322,76 }, { 324,116 }, { 326,156 }, { 328,196 }, { 330,236 }, { 332,276 }, { 334,316 }, { 336,356 }, { 338,396 },
          { 360,38 }, { 362,78 }, { 364,118 }, { 366,158 }, { 368,198 }, { 370,238 }, { 372,278 }, { 374,318 }, { 376,358 }, { 378,398 }
        };
        
        //int syoutotu[1][2] = { {0,20} };
        int KAUNT = 0;
        int j = 0;
        lbfgsfloatval_t R = (lbfgsfloatval_t)RADIUS;
        lbfgsfloatval_t ks = (lbfgsfloatval_t)Ks;
        lbfgsfloatval_t kc = (lbfgsfloatval_t)Kc;

        

        for (int i = 0; i < YARN; i++) {
            t0 = 0.0;
            KAUNT = SEGMENT * 2 * i;
            
            g[SEGMENT * 2 * i] = (ks / 2.0)*((12.0 / 5.0) * x[SEGMENT * 2 * i] + (1.0 / 5.0) * x[SEGMENT * 2 * i + 1] - (12.0 / 5.0) * x[SEGMENT * 2 * i + 2] + (1.0 / 5.0) * x[SEGMENT * 2 * i + 3]);
            g[SEGMENT * 2 * i + 1] = (ks / 2.0)*((1.0 / 5.0) * x[SEGMENT * 2 * i] + (4.0 / 15.0) * x[SEGMENT * 2 * i + 1] - (1.0 / 5.0) * x[SEGMENT * 2 * i + 2] - (1.0 / 15.0) * x[SEGMENT * 2 * i + 3]);
            g[SEGMENT * 2 - 2 + KAUNT] = (ks / 2.0) * (-(12.0 / 5.0) * x[SEGMENT * 2 - 4 + KAUNT] - (1.0 / 5.0) * x[SEGMENT * 2 - 3 + KAUNT] + (12.0 / 5.0) * x[SEGMENT * 2 - 2 + KAUNT] - (1.0 / 5.0) * x[SEGMENT * 2 - 1 + KAUNT]);
            g[SEGMENT * 2 - 1 + KAUNT] = (ks / 2.0) * ((1.0 / 5.0) * x[SEGMENT * 2 - 4 + KAUNT] - (1.0 / 15.0) * x[SEGMENT * 2 - 3 + KAUNT] - (1.0 / 5.0) * x[SEGMENT * 2 - 2 + KAUNT] + (4.0 / 15.0) * x[SEGMENT * 2 - 1 + KAUNT]);
            //printf("%d  %d  \n", SEGMENT * 2 * i, SEGMENT * 2 * i + 1);
            for (j = KAUNT + 2; j < SEGMENT * 2 - 2 + KAUNT; j += 2) {
                g[j] = (ks / 2.0) * (-(12.0 / 5.0) * x[j - 2] - (1.0 / 5.0) * x[j - 1] + (24.0 / 5.0) * x[j] - (12.0 / 5.0) * x[j + 2] + (1.0 / 5.0) * x[j + 3]);
                g[j + 1] = (ks / 2.0) * ((1.0 / 5.0) * x[j - 2] - (1.0 / 15.0) * x[j - 1] + (8.0 / 15.0) * x[j + 1] - (1.0 / 5.0) * x[j + 2] - (1.0 / 15.0) * x[j + 3]);
                //printf("%d  %d\n", j, j + 1);
                t0 += (ks / 2.0) * ((6.0 / 5.0) * x[j] * x[j] + (2.0 / 15.0) * x[j + 1] * x[j + 1]
                    + (6.0 / 5.0) * x[j + 2] * x[j + 2] + (2.0 / 15.0) * x[j + 3] * x[j + 3]
                    + (1.0 / 5.0) * x[j] * x[j + 1] - (12.0 / 5.0) * x[j] * x[j + 2] + (1.0 / 5.0) * x[j] * x[j + 3]
                    - (1.0 / 5.0) * x[j + 2] * x[j + 1] - (1.0 / 15.0) * x[j + 1] * x[j + 3] - (1.0 / 5.0) * x[j + 2] * x[j + 3]);
            }
            //printf("%d  %d  \n", SEGMENT * 2 - 2 + KAUNT, SEGMENT * 2 - 1 + KAUNT);
            t2 += t0;
            //printf("%f %f\n", t0, t2);
            
        }
        for (int l = 0; l < YARN * YARN / 4; l++) {
            lbfgsfloatval_t dif = x[syoutotu[l][0]] - x[syoutotu[l][1]];
            //printf("(%d)dif %f\n", l, dif);
            d = (dif / (R + R)) * (dif / (R + R));
            //printf("%f %f\n", d,dif);
            if (d < 1.0) {
                //printf("衝突起こったよ\n");
                //printf("AAAAAAA\ng[%d] = %f g[%d] = %f\n", syoutotu[l][0], g[syoutotu[l][0]], syoutotu[l][1], g[syoutotu[l][1]]);
                g[syoutotu[l][0]] += kc * (((2.0 * dif) / ((R + R) * (R + R))) - ((2.0 * (R + R) * (R + R)) / (dif * dif * dif)));
                g[syoutotu[l][1]] += kc * (((2.0 * (R + R) * (R + R)) / (dif * dif * dif)) - ((2.0 * dif) / ((R + R) * (R + R))));
                //printf("+ %f + %f\n", kc * (((2.0 * dif) / ((R + R) * (R + R))) - ((2.0 * (R + R) * (R + R)) / (dif * dif * dif))), kc * (((2.0 * (R + R) * (R + R)) / (dif * dif * dif)) - ((2.0 * dif) / ((R + R) * (R + R)))));
                //printf("BBBBBBB\n%d g[%d] = %f g[%d] = %f\n\n",l, syoutotu[l][0], g[syoutotu[l][0]], syoutotu[l][1], g[syoutotu[l][1]]);
                t1 += kc * (1.0 / (d * d) + (d * d) - 2.0);
            }
            else {
                //printf("a");
                t1 += 0.0;
            }

        }
        fx = t1 + t2;
        //printf("%f %f %f\n", t0, t1, fx);

        for (int kkk = 0; kkk < N; kkk = kkk + 2) {
           //printf("x[%d] = %f x[%d] = %f g[%d] = %f g[%d] = %f fx = %f\n", kkk, x[kkk], kkk + 1, x[kkk + 1], kkk, g[kkk], kkk + 1, g[kkk + 1], fx);
        }
        
        for (int l = 0; l < YARN * YARN / 4; l++) {
            //printf("(%d) %d  %d  %f  %f  {%f}  %f  %f\n", l, syoutotu[l][0], syoutotu[l][1], x[syoutotu[l][0]], x[syoutotu[l][1]],x[syoutotu[l][0]] - x[syoutotu[l][1]], g[syoutotu[l][0]], g[syoutotu[l][1]]);
            //printf("%d %d %d\n",l, syoutotu[l][0], syoutotu[l][1]);
        }
        
        return fx;
    }

    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        return reinterpret_cast<objective_function*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        //printf("Iteration %d:\n", k);
        //printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
        //printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
        //printf("\n");
        return 0;
    }
};

lbfgsfloatval_t C0(lbfgsfloatval_t s) {
    return 2 * s * s * s - 3 * s * s + 1;
}

lbfgsfloatval_t C1(lbfgsfloatval_t s) {
    return s * s * s - 2 * s * s + s;
}

lbfgsfloatval_t C2(lbfgsfloatval_t s) {
    return -2 * s * s * s + 3 * s * s;
}

lbfgsfloatval_t C3(lbfgsfloatval_t s) {
    return s * s * s - s * s;
}

void GAISEKI(GLfloat v1[], GLfloat v2[], GLfloat v3[], GLdouble v4[]) {
    GLfloat V1[3];
    GLfloat V2[3];
    for (int i = 0; i < 3; i++) {
        V1[i] = v2[i] - v1[i];
    }
    for (int i = 0; i < 3; i++) {
        V2[i] = v3[i] - v1[i];
    }
    v4[0] = V1[1] * V2[2] - V1[2] * V2[1];
    v4[1] = V1[2] * V2[0] - V1[0] * V2[2];
    v4[2] = V1[0] * V2[1] - V1[1] * V2[0];
}




GLfloat supura[SPLINE][3];
GLfloat supura2[CIRCLE][SPLINE][3];
GLfloat supura2o1[10][SPLINE][3];
GLfloat supura2o2[10][20][SPLINE][3];
GLfloat supura21[CIRCLE][SPLINE][3];
GLfloat supura22[CIRCLE][SPLINE][3];
GLfloat supura3[PLY][SPLINE][3];
GLfloat supura4[PLY][CIRCLE][SPLINE][3];
GLfloat supura4o1[PLY][10][SPLINE][3];
GLfloat supura4o2[PLY][10][20][SPLINE][3];
GLfloat supura41[PLY][CIRCLE][SPLINE][3];
GLfloat supura42[PLY][CIRCLE][SPLINE][3];
GLfloat supura5[PLY][FIBER][SPLINE][3];
GLfloat supura6[PLY][FIBER][CIRCLEF][SPLINE][3];
GLfloat supura61[PLY][FIBER][CIRCLEF][SPLINE][3];
GLfloat supura62[PLY][FIBER][CIRCLEF][SPLINE][3];

GLfloat white[] = { 1.0,1.0,1.0,1.0 };
GLfloat black[] = { 0.0,0.0,0.0,1.0 };
GLfloat half[] = { 0.5,0.5,0.5,1.0 };
GLfloat half2[] = { 0.1,0.1,0.1,1.0 };
GLfloat red[] = { 1.0,0.0,0.0,1.0 };
GLfloat red2[PLY][4];
GLfloat red3[FIBER][4];
GLfloat blue[] = { 0.0,0.0,1.0,1.0 };
GLfloat blue2[PLY][4];
GLfloat blue3[FIBER][4];
GLfloat green[] = { 0.0,1.0,0.0,1.0 };
GLfloat pink[] = { 1.0,0.501,0.501,1.0 };

GLfloat light0pos[] = { 0.0,3.0,-10,0,1.0 };
GLfloat light1pos[] = { 0.0,3.0,-7,0,1.0 };

static GLfloat lightPosition[4] = { 10.0,10.0,7.0,1.0 };
static GLfloat lightDiffuse[3] = {1.0,1.0,1.0 };
static GLfloat lightAmbient[3] = { 0.25,0.25,0.25 };
static GLfloat lightSpecular[3] = { 0.1,0.1,0.1 };
GLfloat spotDirrection[] = { -0.5,0.0,-1.0 };

GLdouble V1[CIRCLE][SPLINE][3];
GLdouble V2[CIRCLE][SPLINE][3];
GLdouble V3[PLY][CIRCLE][SPLINE][3];
GLdouble V4[PLY][CIRCLE][SPLINE][3];
GLdouble V5[PLY][FIBER][CIRCLEF][SPLINE][3];
GLdouble V6[PLY][FIBER][CIRCLEF][SPLINE][3];

void NORMAL2() {
    int last = 0;
    int count = 0;
    GLfloat vv = 0;
    for (int m = 0; m < YARN; m++) {
        for (int k = count; k < SEGMENT + (SEGMENT - 1) * (DIMENSION - 1) - 1 + count; k++) {
            for (int l = 0; l < CIRCLE - 1; l++) {
                
                GAISEKI(&supura2[l][k][0], &supura2[l + 1][k][0], &supura2[l][k + 1][0], &V1[l][k][0]);
                last = l + 1;
                
            }
            GAISEKI(&supura2[last][k][0], &supura2[0][k][0], &supura2[last][k + 1][0], &V1[last][k][0]);
           
        }
        count += SEGMENT * DIMENSION;
    }

        for (int k = 0; k < SPLINE; k++) {
            for (int l = 0; l < CIRCLE; l++) {
                vv = sqrt(V1[l][k][0] * V1[l][k][0] + V1[l][k][1] * V1[l][k][1] + V1[l][k][2] * V1[l][k][2]);
                V2[l][k][0] = V1[l][k][0] / vv;
                V2[l][k][1] = V1[l][k][1] / vv;
                V2[l][k][2] = V1[l][k][2] / vv;
                
            }
        }
    

}

void NORMAL4() {
    int last = 0;
    GLfloat vv = 0;
    for (int i = 0; i < PLY; i++) {
        for (int k = 0; k < SPLINE; k++) {
            for (int l = 0; l < CIRCLE - 1; l++) {
                GAISEKI(&supura4[i][l][k][0], &supura4[i][l + 1][k][0], &supura4[i][l][k + 1][0], &V3[i][l][k][0]);
                last = l + 1;
            }
            GAISEKI(&supura4[i][last][k][0], &supura4[i][0][k][0], &supura4[i][last][k + 1][0], &V3[i][last][k][0]);
        }

    }
    for (int i = 0; i < PLY; i++) {
        for (int k = 0; k < SPLINE; k++) {
            for (int l = 0; l < CIRCLE; l++) {
                vv = sqrt(V3[i][l][k][0] * V3[i][l][k][0] + V3[i][l][k][1] * V3[i][l][k][1] + V3[i][l][k][2] * V3[i][l][k][2]);
                V4[i][l][k][0] = V3[i][l][k][0] / vv;
                V4[i][l][k][1] = V3[i][l][k][1] / vv;
                V4[i][l][k][2] = V3[i][l][k][2] / vv;
            }
        }

    }
}

void NORMAL6() {
    int last = 0;
    GLfloat vv = 0;
    for (int i = 0; i < PLY; i++) {
        for (int j = 0; j < FIBER; j++) {
            for (int k = 0; k < SPLINE; k++) {
                for (int l = 0; l < CIRCLEF - 1; l++) {
                    GAISEKI(&supura6[i][j][l][k][0], &supura6[i][j][l + 1][k][0], &supura6[i][j][l][k + 1][0], &V5[i][j][l][k][0]);
                    last = l + 1;
                }
                GAISEKI(&supura6[i][j][last][k][0], &supura6[i][j][0][k][0], &supura6[i][j][last][k + 1][0], &V5[i][j][last][k][0]);
            }
        }
    }
    for (int i = 0; i < PLY; i++) {
        for (int j = 0; j < FIBER; j++) {
            for (int k = 0; k < SPLINE; k++) {
                for (int l = 0; l < CIRCLEF; l++) {
                    vv = sqrt(V5[i][j][l][k][0] * V5[i][j][l][k][0] + V5[i][j][l][k][1] * V5[i][j][l][k][1] + V5[i][j][l][k][2] * V5[i][j][l][k][2]);
                    V6[i][j][l][k][0] = V5[i][j][l][k][0] / vv;
                    V6[i][j][l][k][1] = V5[i][j][l][k][1] / vv;
                    V6[i][j][l][k][2] = V5[i][j][l][k][2] / vv;
                }
            }
        }
    }
}


void idle(void){
    glutPostRedisplay();
}

void display(void) {
    static double r = 0; 

    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glColor4fv(&white[0]); 
    glLoadIdentity();


    
    gluLookAt(0.0, 0.0, 10.0*RATIO, 0.0, 0.0, 0.0, 0.0, 7.0, 0.0);

    glLightfv(GL_LIGHT0, GL_POSITION, light0pos);
    glLightfv(GL_LIGHT1, GL_POSITION, light1pos);

    glTranslatef(cos(2.0 * r / 360.0 * pi) * -4.5, -4.5,sin(2.0 * r / 360.0 * pi)*4.5);

    glRotated( r, 0.0, 1.0, 0.0);

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, pink);
    /*ヤーンの中心線
    glBegin(GL_LINES);
    int count = 0;
    for (int di = 0; di < YARN; di++) {
        if (di % 2 == 0) {
            glColor4fv(&red[0]);
        }
        else {
            glColor4fv(&blue[0]);
        }
        for (int dj = count; dj < SEGMENT+(SEGMENT - 1)*(DIMENSION-1)-1 + count; dj++) {
            glVertex3fv(&supura[dj][0]);
            glVertex3fv(&supura[dj + 1][0]);
        }
        count += SEGMENT * DIMENSION;
    }
    glEnd();

    */




    /*ヤーンの中心線半径あり
    glBegin(GL_QUADS);
    int count3 = 0;
    for (int di = 0; di < YARN; di++) {
        if (di % 2 == 0) {
            glColor4fv(&red[0]);
            for (int dj = count3; dj < SEGMENT + (SEGMENT - 1) * (DIMENSION - 1) - 1 + count3; dj++) {
                int last = 0;
                for (int dk = 0; dk < circle - 1; dk++) {
                    glVertex3fv(&supura2[dk][dj][0]);
                    glVertex3fv(&supura2[dk][dj + 1][0]);
                    glVertex3fv(&supura2[dk + 1][dj + 1][0]);
                    glVertex3fv(&supura2[dk + 1][dj][0]);
                    last = dk + 1;
                }
                glVertex3fv(&supura2[last][dj][0]);
                glVertex3fv(&supura2[last][dj + 1][0]);
                glVertex3fv(&supura2[0][dj + 1][0]);
                glVertex3fv(&supura2[0][dj][0]);
            }
            count3 += SEGMENT * DIMENSION;
        }
        else {
            glColor4fv(&blue[0]);
            for (int dj = count3; dj < SEGMENT + (SEGMENT - 1) * (DIMENSION - 1) - 1 + count3; dj++) {
                int last = 0;
                for (int dk = 0; dk < circle - 1; dk++) {
                    glVertex3fv(&supura2[dk][dj][0]);
                    glVertex3fv(&supura2[dk + 1][dj][0]);
                    glVertex3fv(&supura2[dk + 1][dj + 1][0]);
                    glVertex3fv(&supura2[dk][dj + 1][0]);
                    last = dk + 1;
                }
                glVertex3fv(&supura2[last][dj][0]);
                glVertex3fv(&supura2[0][dj][0]);
                glVertex3fv(&supura2[0][dj + 1][0]);
                glVertex3fv(&supura2[last][dj + 1][0]);
            }
            count3 += SEGMENT * DIMENSION;
        }
    }
    glEnd();
   // */
   /*ヤーンの全セグメント
   glPointSize(0.01);
   glBegin(GL_POINTS);
   glColor4fv(&white[0]);
   int kaunt = 0;
   for (int dk = 0; dk < YARN; dk++) {
       for (int dl = kaunt; dl < SEGMENT + (SEGMENT - 1) * (DIMENSION - 1) + kaunt; dl++) {
           glVertex3fv(&supura[dl][0]);
       }
       kaunt += SEGMENT * DIMENSION;
   }
   glEnd();
   */
   /*PLYの中心線
   glBegin(GL_LINES);
   int count = 0;
   for (int di = 0; di < YARN; di++) {
       if (di % 2 == 0) {
           glColor4fv(&red[0]);
       }
       else {
           glColor4fv(&blue[0]);
       }
       for (int dk = 0; dk < PLY; dk++) {
           for (int dj = count; dj < SEGMENT + (SEGMENT - 1) * (DIMENSION - 1) - 1 + count; dj++) {
               glVertex3fv(&supura3[dk][dj][0]);
               glVertex3fv(&supura3[dk][dj + 1][0]);
           }
       }
       count += SEGMENT * DIMENSION;
   }
   glEnd();
   // */
   /*PLY半径ありの線(面なし)
   glBegin(GL_LINES);
   int count = 0;
   for (int di = 0; di < YARN; di++) {
       if (di % 2 == 0) {
           glColor4fv(&red[0]);
       }
       else {
           glColor4fv(&blue[0]);
       }
       for (int dh = 0; dh < PLY; dh++) {
           for (int dk = 0; dk < circle; dk++) {
               for (int dj = count; dj < SEGMENT + (SEGMENT - 1) * (DIMENSION - 1) - 1 + count; dj++) {
                   glVertex3fv(&supura4[dh][dk][dj][0]);
                   glVertex3fv(&supura4[dh][dk][dj + 1][0]);
               }
           }
       }
       count += SEGMENT * DIMENSION;
   }
   glEnd();
   // */

   /*PLYの中心線半径あり
   glBegin(GL_QUADS);
   int count3 = 0;
   for (int di = 0; di < YARN; di++) {
       if (di % 2 == 0) {
           for (int dh = 0; dh < PLY; dh++) {
               for (int dc = 0; dc < PLY; dc++) {
                   if (dh == dc) {
                       glColor4fv(&red2[dc][0]);
                   }
               }
               for (int dj = count3; dj < SEGMENT + (SEGMENT - 1) * (DIMENSION - 1) - 1 + count3; dj++) {
                   int last = 0;
                   for (int dk = 0; dk < circle - 1; dk++) {
                       glVertex3fv(&supura4[dh][dk][dj][0]);
                       glVertex3fv(&supura4[dh][dk][dj + 1][0]);
                       glVertex3fv(&supura4[dh][dk + 1][dj + 1][0]);
                       glVertex3fv(&supura4[dh][dk + 1][dj][0]);
                       last = dk + 1;
                   }
                   glVertex3fv(&supura4[dh][last][dj][0]);
                   glVertex3fv(&supura4[dh][last][dj + 1][0]);
                   glVertex3fv(&supura4[dh][0][dj + 1][0]);
                   glVertex3fv(&supura4[dh][0][dj][0]);
               }

           }
           count3 += SEGMENT * DIMENSION;
       }
       else {

           for (int dh = 0; dh < PLY; dh++) {
               for (int dc = 0; dc < PLY; dc++) {
                   if (dh == dc) {
                       glColor4fv(&blue2[dc][0]);
                   }
               }
               for (int dj = count3; dj < SEGMENT + (SEGMENT - 1) * (DIMENSION - 1) - 1 + count3; dj++) {
                   int last = 0;
                   for (int dk = 0; dk < circle - 1; dk++) {
                       glVertex3fv(&supura4[dh][dk][dj][0]);
                       glVertex3fv(&supura4[dh][dk + 1][dj][0]);
                       glVertex3fv(&supura4[dh][dk + 1][dj + 1][0]);
                       glVertex3fv(&supura4[dh][dk][dj + 1][0]);
                       last = dk + 1;
                   }
                   glVertex3fv(&supura4[dh][last][dj][0]);
                   glVertex3fv(&supura4[dh][0][dj][0]);
                   glVertex3fv(&supura4[dh][0][dj + 1][0]);
                   glVertex3fv(&supura4[dh][last][dj + 1][0]);
               }

           }
           count3 += SEGMENT * DIMENSION;
       }

   }
   glEnd();
    //*/

    /*FIBERの中心
    glBegin(GL_LINES);
    int count = 0;
    for (int di = 0; di < YARN; di++) {
        if (di % 2 == 0) {
            glColor4fv(&red[0]);
        }
        else {
            glColor4fv(&blue[0]);
        }
        for (int dk = 0; dk < PLY; dk++) {
            for (int df = 0; df < FIBER; df++) {
                for (int dj = count; dj < SEGMENT + (SEGMENT - 1) * (DIMENSION - 1) - 1 + count; dj++) {
                    glVertex3fv(&supura5[dk][df][dj][0]);
                    glVertex3fv(&supura5[dk][df][dj + 1][0]);
                }
            }
        }
        count += SEGMENT * DIMENSION;
    }
    glEnd();
    // */

    /*FIBER半径ありの線(面なし)
    glBegin(GL_LINES);
    int count = 0;
    for (int di = 0; di < YARN; di++) {
        if (di % 2 == 0) {
            glColor4fv(&red[0]);
        }
        else {
            glColor4fv(&blue[0]);
        }
        for (int dh = 0; dh < PLY; dh++) {
            for (int df = 0; df < FIBER; df++) {
                for (int dk = 0; dk < circle; dk++) {
                    for (int dj = count; dj < SEGMENT + (SEGMENT - 1) * (DIMENSION - 1) - 1 + count; dj++) {
                        glVertex3fv(&supura6[dh][df][dk][dj][0]);
                        glVertex3fv(&supura6[dh][df][dk][dj + 1][0]);
                    }
                }
            }
        }
        count += SEGMENT * DIMENSION;
    }
    glEnd();
    // */

    /*糸光点あり
    glBegin(GL_QUADS);
    int count3 = 0;
    for (int di = 0; di < YARN; di++) {
        for (int dj = count3; dj < (SEGMENT - 1) * DIMENSION -1  + count3; dj++) {
            int last = 0;
            for (int dk = 0; dk < CIRCLE - 1; dk++) {
                glNormal3dv(&V2[dk][dj][0]);
                glVertex3fv(&supura2[dk][dj][0]);
                glVertex3fv(&supura2[dk + 1][dj][0]);
                glVertex3fv(&supura2[dk + 1][dj + 1][0]);
                glVertex3fv(&supura2[dk][dj + 1][0]);
                last = dk + 1;
            }
            
            glNormal3dv(&V2[last][dj][0]);
            glVertex3fv(&supura2[last][dj][0]);
            glVertex3fv(&supura2[0][dj][0]);
            glVertex3fv(&supura2[0][dj + 1][0]);
            glVertex3fv(&supura2[last][dj + 1][0]);
        }
        count3 += SEGMENT * DIMENSION;


    }
    glEnd();
    //*/

/*PLY光点
   glBegin(GL_QUADS);
   int count3 = 0;
   for (int di = 0; di < YARN; di++) {

       for (int dh = 0; dh < PLY; dh++) {

           for (int dj = count3; dj < SEGMENT + (SEGMENT - 1) * (DIMENSION - 1) - 1 + count3; dj++) {
               int last = 0;
               for (int dk = 0; dk < circle - 1; dk++) {
                   glNormal3dv(&V4[dh][dk][dj][0]);
                   glVertex3fv(&supura4[dh][dk][dj][0]);
                   glVertex3fv(&supura4[dh][dk + 1][dj][0]);
                   glVertex3fv(&supura4[dh][dk + 1][dj + 1][0]);
                   glVertex3fv(&supura4[dh][dk][dj + 1][0]);
                   last = dk + 1;
               }
               
               glNormal3dv(&V4[dh][last][dj][0]);
               glVertex3fv(&supura4[dh][last][dj][0]);
               glVertex3fv(&supura4[dh][0][dj][0]);
               glVertex3fv(&supura4[dh][0][dj + 1][0]);
               glVertex3fv(&supura4[dh][last][dj + 1][0]);
           }

       }
       count3 += SEGMENT * DIMENSION;



   }
   glEnd();
    //*/


    //*FIBERの中心線半径あり
    glBegin(GL_QUADS);
    int count3 = 0;
    for (int di = 0; di < YARN; di++) {

        for (int dh = 0; dh < PLY; dh++) {
            for (int df = 0; df < FIBER; df++) {
                for (int dj = count3; dj < SEGMENT + (SEGMENT - 1) * (DIMENSION - 1) - 1 + count3; dj++) {
                    int last = 0;
                    for (int dk = 0; dk < CIRCLEF - 1; dk++) {
                       glNormal3dv(&V6[dh][df][dk][dj][0]);
                        glVertex3fv(&supura6[dh][df][dk][dj][0]);
                        glVertex3fv(&supura6[dh][df][dk + 1][dj][0]);
                        glVertex3fv(&supura6[dh][df][dk + 1][dj + 1][0]);
                        glVertex3fv(&supura6[dh][df][dk][dj + 1][0]);
                        last = dk + 1;
                    }
                    glNormal3dv(&V6[dh][df][last][dj][0]);
                    glVertex3fv(&supura6[dh][df][last][dj][0]);
                    glVertex3fv(&supura6[dh][df][0][dj][0]);
                    glVertex3fv(&supura6[dh][df][0][dj + 1][0]);
                    glVertex3fv(&supura6[dh][df][last][dj + 1][0]);
                }
            }
        }
        count3 += SEGMENT * DIMENSION;

    }
    glEnd();
    //*/
    
    /*
    glPointSize(0.01);
    glBegin(GL_POINTS);
    glColor4fv(&white[0]);
    int kaunt = 0;
    for (int dk = 0; dk < YARN; dk++) {
        for (int dl = kaunt; dl < SEGMENT + (SEGMENT - 1) * (DIMENSION - 1) + kaunt; dl++) {
            glVertex3fv(&supura[dl][0]);
        }
        kaunt += SEGMENT * DIMENSION;
    }
    glEnd();
    //*/
    /*
    glBegin(GL_POINTS);
    glColor4fv(&red[0]);
    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 10; j++) {
            glVertex3f(supura6[1][j][i][1000][0],supura6[1][j][i][1000][1],supura6[1][j][i][1000][2]);
        }
    }
    glEnd();
    //*/
/*
    glBegin(GL_POINTS);
    glColor4fv(&white[0]);
    for (int i = 0; i < 100; i++) {
            glVertex3f(supura4[1][i][800][0],supura4[1][i][800][1], supura4[1][i][800][2]);
        
    }
    glEnd();
    /*
    glBegin(GL_POINTS);
    glColor4fv(&white[0]);
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            glVertex2f(O[i][j][0], O[i][j][1]);
        }
    }
    glEnd();
    //*/

    glutIdleFunc(idle);
    

    glutSwapBuffers();
    r += 0.1;
    if (r >= 360) r = 0;
}

void resize(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(100.0, (double)w / (double)h, 1.0, 100.0); 
    glMatrixMode(GL_MODELVIEW);
}

void mouse(int button, int state, int x, int y)
{
    switch (button) {
    case GLUT_LEFT_BUTTON:
        /* アニメーション開始 */
        if (state == GLUT_DOWN) {
            glutIdleFunc(idle);
        }
        else {
            glutIdleFunc(0);
        }
        break;
    case GLUT_RIGHT_BUTTON:
        if (state == GLUT_DOWN) {
            /* コマ送り (1ステップだけ進める) */
            glutPostRedisplay();
        }
        break;
    default:
        break;
    }
}


void init(void)
{
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glEnable(GL_DEPTH_TEST);

    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);

    glEnable(GL_LIGHTING);
    //glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    //glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, lightDiffuse);
    glLightfv(GL_LIGHT1, GL_SPECULAR, lightSpecular);
    glLightfv(GL_LIGHT1, GL_AMBIENT, lightAmbient);
    //glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.01);
    //glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, spotDirrection);
    //glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 30.0);
    //glLightf(GL_LIGHT0, GL_SPOT_EXPONENT, 0.0);

    //glLightfv(GL_LIGHT1, GL_DIFFUSE, white);
    //glLightfv(GL_LIGHT1, GL_SPECULAR, half);
    //glLightfv(GL_LIGHT1, GL_AMBIENT, black);
}

void Y_y(int k,lbfgsfloatval_t s) {
    for (int p = 0; p < 10; p++) {
        supura2o1[p][k][0] = supura[k][0];
        supura2o1[p][k][1] = supura[k][1] + radius*sin(p * pi / 5) * 0.75;
        supura2o1[p][k][2] = supura[k][2] + radius*cos(p * pi / 5) * 0.75;
        for (int j = 0; j < 20; j++) {
            supura2o2[p][j][k][0] = supura2o1[p][k][0];
            supura2o2[p][j][k][1] = supura2o1[p][k][1] + radius * 0.25 * sin(j * pi / 10 - (5 - 2 * p) * pi / 10);
            supura2o2[p][j][k][2] = supura2o1[p][k][2] + radius * 0.25 * cos(j * pi / 10 - (5 - 2 * p) * pi / 10);
        }
    }

    for (int i = 0; i < CIRCLE; i++) {
        if (z == 10) {
            z = 0;
            w++;
        }
        if (w == 10) {
            w = 0;
        }
        
        supura2[i][k][0] = supura2o2[z][w][k][0];
        supura2[i][k][1] = supura2o2[z][w][k][1];
        supura2[i][k][2] = supura2o2[z][w][k][2];
        z++;
    }
    
    lbfgsfloatval_t r = 0;
    if (k % (DIMENSION*SEGMENT) == 0) {
        r = 0;
    }
    else {
        r = atan2(s, supura[k][2] - supura[k - 1][2]);
        for (int y = 0; y < circle; y++) {
            supura21[y][k][0] = supura2[y][k][0] - supura[k][0];
            supura21[y][k][1] = supura2[y][k][1] - supura[k][1];
            supura21[y][k][2] = supura2[y][k][2] - supura[k][2];

            supura22[y][k][0] = supura21[y][k][0] * cos(pi / 2 - r) - supura21[y][k][2] * sin(pi / 2 - r);
            supura22[y][k][1] = supura21[y][k][1];
            supura22[y][k][2] = supura21[y][k][0] * sin(pi / 2 - r) + supura21[y][k][2] * cos(pi / 2 - r);

            supura2[y][k][0] = supura22[y][k][0] + supura[k][0];
            supura2[y][k][1] = supura22[y][k][1] + supura[k][1];
            supura2[y][k][2] = supura22[y][k][2] + supura[k][2];
        }
    }
}

void Y_t(int k,lbfgsfloatval_t s) {
    for (int p = 0; p < 10; p++) {
        supura2o1[p][k][0] = supura[k][0] + radius * cos(p * pi / 5) * 0.75;
        supura2o1[p][k][1] = supura[k][1];
        supura2o1[p][k][2] = supura[k][2] + radius * sin(p * pi / 5) * 0.75;
        for (int j = 0; j < 20; j++) {
            supura2o2[p][j][k][0] = supura2o1[p][k][0] + radius * 0.25 * cos(j * pi / 10 - (5 - 2 * p) * pi / 10);
            supura2o2[p][j][k][1] = supura2o1[p][k][1];
            supura2o2[p][j][k][2] = supura2o1[p][k][2] + radius * 0.25 * sin(j * pi / 10 - (5 - 2 * p) * pi / 10);
        }
    }
    for (int i = 0; i < CIRCLE; i++) {
        if (z == 10) {
            z = 0;
            w++;
        }
        if (w == 10) {
            w = 0;
        }

        supura2[i][k][0] = supura2o2[z][w][k][0];
        supura2[i][k][1] = supura2o2[z][w][k][1];
        supura2[i][k][2] = supura2o2[z][w][k][2];
        z++;
    
    }
    lbfgsfloatval_t r = 0;
    if (k % (DIMENSION*SEGMENT) == 0) {
        r = 0;
    }
    else {
        r = atan2(s, supura[k][2] - supura[k - 1][2]);
        for (int y = 0; y < circle; y++) {
            supura21[y][k][0] = supura2[y][k][0] - supura[k][0];
            supura21[y][k][1] = supura2[y][k][1] - supura[k][1];
            supura21[y][k][2] = supura2[y][k][2] - supura[k][2];

            supura22[y][k][0] = supura21[y][k][0];
            supura22[y][k][1] = supura21[y][k][1] * cos(pi / 2 - r) - supura21[y][k][2] * sin(pi / 2 - r);
            supura22[y][k][2] = supura21[y][k][1] * sin(pi / 2 - r) + supura21[y][k][2] * cos(pi / 2 - r);

            supura2[y][k][0] = supura22[y][k][0] + supura[k][0];
            supura2[y][k][1] = supura22[y][k][1] + supura[k][1];
            supura2[y][k][2] = supura22[y][k][2] + supura[k][2];
        }
    }
}


void P_y(int k) {
    for (int w = 0; w < PLY; w++) {
        lbfgsfloatval_t Iply = Tply + 2 * pi * (w - 1) / PLY;
        lbfgsfloatval_t Sx = Rply * cos(2 * pi * supura[k][0] / Aply + Iply);
        lbfgsfloatval_t Sy = Rply * sin(2 * pi * supura[k][0] / Aply + Iply);
        supura3[w][k][0] = supura[k][0];
        supura3[w][k][1] = supura[k][1] + Sy;
        supura3[w][k][2] = supura[k][2] + Sx;
    }
}

void P_t(int k) {
    for (int w = 0; w < PLY; w++) {
        lbfgsfloatval_t Iply = Tply + 2 * pi * (w - 1) / PLY;
        lbfgsfloatval_t Sx = Rply * cos(2 * pi * supura[k][1] / Aply + Iply);
        lbfgsfloatval_t Sy = Rply * sin(2 * pi * supura[k][1] / Aply + Iply);
        supura3[w][k][0] = supura[k][0] + Sx;
        supura3[w][k][1] = supura[k][1];
        supura3[w][k][2] = supura[k][2] + Sy;
    }
}

void Q_y(int k,lbfgsfloatval_t s) {
    for (int g = 0; g < PLY; g++) {
        for (int h = 0; h < 10; h++) {
            supura4o1[g][h][k][0] = supura3[g][k][0];
            supura4o1[g][h][k][1] = supura3[g][k][1] + YPdistance * sin(k  * pi/20 + h * pi / 5) * 0.75;
            supura4o1[g][h][k][2] = supura3[g][k][2] + YPdistance * cos(k  * pi/20 + h * pi / 5) * 0.75;
            for (int j = 0; j < 20; j++) {
                supura4o2[g][h][j][k][0] = supura4o1[g][h][k][0];
                supura4o2[g][h][j][k][1] = supura4o1[g][h][k][1] + YPdistance * 0.3 * sin(j * pi / 10 - (5 - 2 * h) * pi / 10);
                supura4o2[g][h][j][k][2] = supura4o1[g][h][k][2] + YPdistance * 0.3 * cos(j * pi / 10 - (5 - 2 * h) * pi / 10);
            }
        }
    }
    for (int g = 0; g < PLY; g++) {
        z = 0;
        w = 0;
        for (int i = 0; i < CIRCLE; i++) {
            if (z == 10) {
                z = 0;
                w++;
            }
            if (w == 10) {
                w = 0;
            }

            supura4[g][i][k][0] = supura4o2[g][z][w][k][0];
            supura4[g][i][k][1] = supura4o2[g][z][w][k][1];
            supura4[g][i][k][2] = supura4o2[g][z][w][k][2];
            z++;
        }
    }
   
    lbfgsfloatval_t r = 0;
    if (k % (DIMENSION * SEGMENT) == 0) {
        r = 0;
    }
    else {
        r = atan2(s, supura3[0][k][2] - supura3[0][k - 1][2]);
        for (int x = 0; x < PLY; x++) {
            for (int y = 0; y < circle; y++) {
                supura41[x][y][k][0] = supura4[x][y][k][0] - supura3[x][k][0];
                supura41[x][y][k][1] = supura4[x][y][k][1] - supura3[x][k][1];
                supura41[x][y][k][2] = supura4[x][y][k][2] - supura3[x][k][2];

                supura42[x][y][k][0] = supura41[x][y][k][0] * cos(pi / 2 - r) - supura41[x][y][k][2] * sin(pi / 2 - r);
                supura42[x][y][k][1] = supura41[x][y][k][1];
                supura42[x][y][k][2] = supura41[x][y][k][0] * sin(pi / 2 - r) + supura41[x][y][k][2] * cos(pi / 2 - r);

                supura4[x][y][k][0] = supura42[x][y][k][0] + supura3[x][k][0];
                supura4[x][y][k][1] = supura42[x][y][k][1] + supura3[x][k][1];
                supura4[x][y][k][2] = supura42[x][y][k][2] + supura3[x][k][2];
            }
        }
    }
}

void Q_t(int k ,lbfgsfloatval_t s) {
    for (int g = 0; g < PLY; g++) {
        for (int h = 0; h < 10; h++) {
            supura4o1[g][h][k][0] = supura3[g][k][0] + YPdistance * cos(k *pi/20 + h * pi / 5) * 0.75;
            supura4o1[g][h][k][1] = supura3[g][k][1];
            supura4o1[g][h][k][2] = supura3[g][k][2] + YPdistance * sin(k *pi/20 + h * pi / 5) * 0.75;
            for (int j = 0; j < 20; j++) {
                supura4o2[g][h][j][k][0] = supura4o1[g][h][k][0]+ YPdistance * 0.3 * cos(j * pi / 10 - (5 - 2 * h) * pi / 10);
                supura4o2[g][h][j][k][1] = supura4o1[g][h][k][1];
                supura4o2[g][h][j][k][2] = supura4o1[g][h][k][2] + YPdistance * 0.3 * sin(j * pi / 10 - (5 - 2 * h) * pi / 10);
            }
        }
    }
    for (int g = 0; g < PLY; g++) {
        z = 0;
        w = 0;
        for (int i = 0; i < CIRCLE; i++) {
            if (z == 10) {
                z = 0;
                w++;
            }
            if (w == 10) {
                w = 0;
            }

            supura4[g][i][k][0] = supura4o2[g][z][w][k][0];
            supura4[g][i][k][1] = supura4o2[g][z][w][k][1];
            supura4[g][i][k][2] = supura4o2[g][z][w][k][2];
            z++;
        }
    }
    lbfgsfloatval_t r = 0;
    if (k % (DIMENSION * SEGMENT) == 0) {
        r = 0;
    }
    else {
        r = atan2(s, supura3[0][k][2] - supura3[0][k - 1][2]);
        for (int x = 0; x < PLY; x++) {
            for (int y = 0; y < circle; y++) {
                supura41[x][y][k][0] = supura4[x][y][k][0] - supura3[x][k][0];
                supura41[x][y][k][1] = supura4[x][y][k][1] - supura3[x][k][1];
                supura41[x][y][k][2] = supura4[x][y][k][2] - supura3[x][k][2];

                supura42[x][y][k][0] = supura41[x][y][k][0];
                supura42[x][y][k][1] = supura41[x][y][k][1] * cos(pi / 2 - r) - supura41[x][y][k][2] * sin(pi / 2 - r);
                supura42[x][y][k][2] = supura41[x][y][k][1] * sin(pi / 2 - r) + supura41[x][y][k][2] * cos(pi / 2 - r);

                supura4[x][y][k][0] = supura42[x][y][k][0] + supura3[x][k][0];
                supura4[x][y][k][1] = supura42[x][y][k][1] + supura3[x][k][1];
                supura4[x][y][k][2] = supura42[x][y][k][2] + supura3[x][k][2];
            }
        }
    }
}

void F_y(int k) {
    for (int g = 0; g < PLY; g++) {
        for (int f = 0; f < FIBER; f++) {
            lbfgsfloatval_t Ric = 0.75 * cos(supura3[g][k][0] * 2 * pi * Afiber);
            lbfgsfloatval_t Ris = 0.75 * sin(supura3[g][k][0] * 2 * pi * Afiber);
            supura5[g][f][k][0] = supura3[g][k][0];
            supura5[g][f][k][1] = supura3[g][k][1] + YPdistance * sin(k*2 * pi * Afiber + f * pi / 5) * 0.75;
            supura5[g][f][k][2] = supura3[g][k][2] + YPdistance * cos(k*2 * pi * Afiber + f * pi / 5) * 0.75;
        }
    }
}

void F_t(int k) {
    for (int g = 0; g < PLY; g++) {
        for (int f = 0; f < FIBER; f++) {
            lbfgsfloatval_t Ric = 0.75 * cos(supura3[g][k][1] * 2 * pi * Afiber);
            lbfgsfloatval_t Ris = 0.75 * sin(supura3[g][k][1] * 2 * pi * Afiber);
            supura5[g][f][k][0] = supura3[g][k][0] + YPdistance * cos(k *2 * pi * Afiber + f * pi / 5) * 0.75;
            supura5[g][f][k][1] = supura3[g][k][1];
            supura5[g][f][k][2] = supura3[g][k][2] + YPdistance * sin(k  * 2 * pi * Afiber + f * pi / 5) * 0.75;
        }
    }
}

void G_y(int k,lbfgsfloatval_t s) {
    for (int g = 0; g < PLY; g++) {
        for (int f = 0; f < FIBER; f++) {
            for (int h = 0; h < CIRCLEF; h++) {
                lbfgsfloatval_t Gsin = sin((2.0 / circlef) * h * pi);
                lbfgsfloatval_t Gcos = cos((2.0 / circlef) * h * pi);
                supura6[g][f][h][k][0] = supura5[g][f][k][0];
                supura6[g][f][h][k][1] = supura5[g][f][k][1] + YPdistance * 0.3 * Gsin;
                supura6[g][f][h][k][2] = supura5[g][f][k][2] + YPdistance * 0.3 * Gcos;
            }
        }
    }
    lbfgsfloatval_t r = 0;
    if (k % (DIMENSION * SEGMENT) == 0) {
        r = 0;
    }
    else {
        r = atan2(s, supura5[0][0][k][2] - supura5[0][0][k - 1][2]);
        for (int x = 0; x < PLY; x++) {
            for (int z = 0; z < FIBER; z++) {
                for (int y = 0; y < CIRCLEF; y++) {
                    supura61[x][z][y][k][0] = supura6[x][z][y][k][0] - supura5[x][z][k][0];
                    supura61[x][z][y][k][1] = supura6[x][z][y][k][1] - supura5[x][z][k][1];
                    supura61[x][z][y][k][2] = supura6[x][z][y][k][2] - supura5[x][z][k][2];

                    supura62[x][z][y][k][0] = supura61[x][z][y][k][0] * cos(pi / 2 - r) - supura61[x][z][y][k][2] * sin(pi / 2 - r);
                    supura62[x][z][y][k][1] = supura61[x][z][y][k][1];
                    supura62[x][z][y][k][2] = supura61[x][z][y][k][0] * sin(pi / 2 - r) + supura61[x][z][y][k][2] * cos(pi / 2 - r);

                    supura6[x][z][y][k][0] = supura62[x][z][y][k][0] + supura5[x][z][k][0];
                    supura6[x][z][y][k][1] = supura62[x][z][y][k][1] + supura5[x][z][k][1];
                    supura6[x][z][y][k][2] = supura62[x][z][y][k][2] + supura5[x][z][k][2];
                }
            }
        }
    }
}

void G_t(int k,lbfgsfloatval_t s) {
    for (int g = 0; g < PLY; g++) {
        for (int f = 0; f < FIBER; f++) {
            for (int h = 0; h < CIRCLEF; h++) {
                lbfgsfloatval_t Gsin = sin((2.0 / circlef) * h * pi);
                lbfgsfloatval_t Gcos = cos((2.0 / circlef) * h * pi);
                supura6[g][f][h][k][0] = supura5[g][f][k][0] + YPdistance * 0.3 * Gcos;
                supura6[g][f][h][k][1] = supura5[g][f][k][1];
                supura6[g][f][h][k][2] = supura5[g][f][k][2] + YPdistance * 0.3 * Gsin;
            }
        }
    }
    lbfgsfloatval_t r = 0;
    if (k % (DIMENSION * SEGMENT) == 0) {
        r = 0;
    }
    else {
        r = atan2(s, supura5[0][0][k][2] - supura5[0][0][k - 1][2]);
        for (int x = 0; x < PLY; x++) {
            for (int z = 0; z < FIBER; z++) {
                for (int y = 0; y < CIRCLEF; y++) {
                    supura61[x][z][y][k][0] = supura6[x][z][y][k][0] - supura5[x][z][k][0];
                    supura61[x][z][y][k][1] = supura6[x][z][y][k][1] - supura5[x][z][k][1];
                    supura61[x][z][y][k][2] = supura6[x][z][y][k][2] - supura5[x][z][k][2];

                    supura62[x][z][y][k][0] = supura61[x][z][y][k][0];
                    supura62[x][z][y][k][1] = supura61[x][z][y][k][1] * cos(pi / 2 - r) - supura61[x][z][y][k][2] * sin(pi / 2 - r);
                    supura62[x][z][y][k][2] = supura61[x][z][y][k][1] * sin(pi / 2 - r) + supura61[x][z][y][k][2] * cos(pi / 2 - r);

                    supura6[x][z][y][k][0] = supura62[x][z][y][k][0] + supura5[x][z][k][0];
                    supura6[x][z][y][k][1] = supura62[x][z][y][k][1] + supura5[x][z][k][1];
                    supura6[x][z][y][k][2] = supura62[x][z][y][k][2] + supura5[x][z][k][2];
                }
            }
        }
    }
}

void SEG_y(lbfgsfloatval_t s, int i, GLfloat supura[][3], lbfgsfloatval_t j, int k,lbfgsfloatval_t r) {
    supura[k][0] = j * RATIO;
    supura[k][1] = r;
    supura[k][2] = K[i];
    
    Y_y(k,s);
    P_y(k);
    Q_y(k,s);
    F_y(k);
    G_y(k,s);
    
}

void SPR_y(lbfgsfloatval_t s, lbfgsfloatval_t S,int i,GLfloat supura[][3],lbfgsfloatval_t j,int k, lbfgsfloatval_t r) {
    j = j + s;
    for (int l = 1; l < DIMENSION; l++) {
        double p = s * l;
        S = C0(p) * K[i] + C1(p) * K[i + 1] + C2(p) * K[i + 2] + C3(p) * K[i + 3];
        supura[k][0] = j * RATIO;
        supura[k][1] = r;
        supura[k][2] = S;
        Y_y(k,s);
        P_y(k);
        Q_y(k,s);
        F_y(k);
        G_y(k,s);
        j = j + s;
        k++;
    }
}

void SEG_t(lbfgsfloatval_t s, int i, GLfloat supura[][3], lbfgsfloatval_t j, int k, lbfgsfloatval_t r) {
    supura[k][0] = r;
    supura[k][1] = j * RATIO;
    supura[k][2] = K[i];
    Y_t(k,s);
    P_t(k);
    Q_t(k,s);
    F_t(k);
    G_t(k,s);
}

void SPR_t(lbfgsfloatval_t s, lbfgsfloatval_t S, int i, GLfloat supura[][3], lbfgsfloatval_t j, int k, lbfgsfloatval_t r) {
    j = j + s;
    for (int l = 1; l < DIMENSION; l++) {
        double p = s * l;
        S = C0(p) * K[i] + C1(p) * K[i + 1] + C2(p) * K[i + 2] + C3(p) * K[i + 3];
        supura[k][0] = r;
        supura[k][1] = j * RATIO;
        supura[k][2] = S;
        Y_t(k,s);
        P_t(k);
        Q_t(k,s);
        F_t(k);
        G_t(k,s);
        j = j + s;
        k++;
    }
}

lbfgsfloatval_t p(lbfgsfloatval_t x, lbfgsfloatval_t E, lbfgsfloatval_t Beta) {
    return ((1.0 - 2.0 * E) * pow((exp(1) - exp(x)) / (exp(1) - 1.0), Beta) + E);
}

lbfgsfloatval_t sampling(lbfgsfloatval_t E, lbfgsfloatval_t Beta, lbfgsfloatval_t* en) {
    lbfgsfloatval_t z1 = 10.0;
    lbfgsfloatval_t z2 = 10.0;
    lbfgsfloatval_t z3 = 10.0;
    lbfgsfloatval_t z4 = 10.0;
    lbfgsfloatval_t u = 0.0;
    std::random_device rnd;     // 非決定的な乱数生成器を生成
    std::mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    std::uniform_real_distribution<> rand100(-1, 1);        // [0, 99] 範囲の一様乱数
    std::uniform_real_distribution<> rand101(0, 1);
    while (1) {
        while (1) {
            z1 = rand100(mt);
            z2 = rand100(mt);
            if (sqrt(z1 * z1 + z2 * z2) < 1.0) {
                break;
            }
        }
        while (1) {
            z3 = rand100(mt);
            z4 = rand100(mt);
            if (sqrt(z3 * z3 + z4 * z4) < 1.0) {
                break;
            }
        }
        E = z1;
        Beta = z2;
        u = rand101(mt);
        //printf("%f %f ", z1, z2);
        //printf("%f %f\n", p(sqrt(z1 * z1 + z2 * z2), E, Beta), u);
        if (p(sqrt(z3 * z3 + z4 * z4), E, Beta) > u) {
            en[0] = z1;
            en[1] = z2;
            return 0;
        }
    }
}

int main(int argc, char **argv)
{
    objective_function obj;
    obj.run(N);
    lbfgsfloatval_t S = 0.0;
    lbfgsfloatval_t s = 1.0 / (lbfgsfloatval_t)DIMENSION;
    lbfgsfloatval_t j = 0.0;
    lbfgsfloatval_t interval = (lbfgsfloatval_t)INTERVAL;
    lbfgsfloatval_t V = 0.0;
    lbfgsfloatval_t E = 0.0;
    lbfgsfloatval_t Beta = 0.0;
    lbfgsfloatval_t EN[2];
    lbfgsfloatval_t* en;
    en = &EN[0];
    int k = 0;
    int COUNT = 0;
    for (int i = 0; i < 10; i++) {
        o[i][0] = 10.0 + cos(i * pi / 5) * 0.75;
        o[i][1] = 10.0  + sin(i * pi / 5) *0.75;
        for (int j = 0; j < 20; j++) {
            O[i][j][0] = o[i][0] + 0.3 * cos(j * pi / 10-(5 - 2*i) * pi / 10);
            O[i][j][1] = o[i][1] + 0.3 * sin(j * pi / 10-(5 - 2*i) * pi / 10);
            
        }
        
    }
    /*
    for (int i = 0; i < FIBER; i++) {
        sampling(E, Beta, en);
        F[i][0] = ;
        F[i][1] = ;
        printf("%f %f\n", EN[0], EN[1]);
    }
    */
    for (int m = 0; m < YARN; m++) {
        if (m % 2 == 0) {
            for (int i = COUNT; i < N / YARN + COUNT; i = i + 2) {
                SEG_y(s,i, supura, j, k, V);
                k++;
                SPR_y(s, S, i, supura, j, k, V);
                j++;
                k = k + DIMENSION - 1;
            }
        }
        else {
            for (int i = COUNT; i < N / YARN + COUNT; i = i + 2) {
                SEG_t(s,i, supura, j, k, V);
                k++;
                SPR_t(s, S, i, supura, j, k, V);
                j++;
                k = k + DIMENSION - 1;
            }
            V += interval;
        }
        COUNT += N / YARN;
        j = 0;
    }

    for (int i = 0; i < SPLINE; i++) {
        for (int j = 0; j < CIRCLE; j++) {
            //printf("%d %d %f %f %f \n",i , j, supura2[i][j][0], supura2[i][j][1], supura2[i][j][2]);
        }
    }

    lbfgsfloatval_t ply = (lbfgsfloatval_t)PLY;
    GLfloat R = 0.0, G = 0.0, B = 0.0;

    R = 1.0;
    for (int c = 0; c < PLY; c++) {
        red2[c][0] = R;
        red2[c][1] = G;
        red2[c][2] = B;
        red2[c][3] = 1.0;
        G += 1.0 / ply;
        B += 1.0 / ply;
    }
    R = 0.0; G = 0.0; B = 0.0;
    B = 1.0;
    for (int c = 0; c < PLY; c++) {
        blue2[c][0] = R;
        blue2[c][1] = G;
        blue2[c][2] = B;
        blue2[c][3] = 1.0;
        G += 1.0 / ply;
        R += 1.0 / ply;
    }

    GLfloat Rf = 0.0, Gf = 0.0, Bf = 0.0;
    Rf = 1.0;
    for (int c = 0; c < FIBER; c++) {
        red3[c][0] = Rf;
        red3[c][1] = Gf;
        red3[c][2] = Bf;
        red3[c][3] = 1.0;
        Gf += 1.0 / (GLfloat)FIBER;
        Bf += 1.0 / (GLfloat)FIBER;
    }
    Rf = 0.0; Gf = 0.0; Bf = 0.0;
    Bf = 1.0;
    for (int c = 0; c < FIBER; c++) {
        blue3[c][0] = Rf;
        blue3[c][1] = Gf;
        blue3[c][2] = Bf;
        blue3[c][3] = 1.0;
        Gf += 1.0 / (GLfloat)FIBER;
        Rf += 1.0 / (GLfloat)FIBER;
    }
    NORMAL2();
    NORMAL4();
    NORMAL6();

    
    glutInit(&argc, argv);
    glutInitWindowPosition(100, 50);
    glutInitWindowSize(700, 700);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("Hello GLUT!!");
    glutDisplayFunc(display);
    glutReshapeFunc(resize);
    glutMouseFunc(mouse);
    init();
    glutMainLoop();
    return 0;
}
