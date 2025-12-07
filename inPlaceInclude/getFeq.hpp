// Calculates floats feq0, feq1 ... from floats rho, ux, uy, uz
// Floats rho, ux, uy, uz must be defined above

// id: { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
// cx: { 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
// cy: { 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
// cz: { 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };
// w:  { 8/27, 2/27, 2/27, 2/27 , 2/27, 2/27, 2/27, 1/54, 1/54, 1/54, 1/54, 1/54, 1/54, 1/54, 1/54, 1/54, 1/54, 1/54, 1/54, 1/216, 1/216, 1/216, 1/216, 1/216, 1/216, 1/216, 1/216 };

const float u2 = ux*ux + uy*uy + uz*uz;

const float cu0  = 0.f;
const float cu1  = +ux;
const float cu2  = -ux;
const float cu3  = -uz;
const float cu4  = +uz;
const float cu5  = -uy;
const float cu6  = +uy;
const float cu7  = +ux -uz;
const float cu8  = -ux +uz;
const float cu9  = +ux +uz;
const float cu10 = -ux -uz;
const float cu11 = -ux -uy;
const float cu12 = +ux +uy;
const float cu13 = +uy -uz;
const float cu14 = -uy +uz;
const float cu15 = -ux +uy;
const float cu16 = +ux -uy;
const float cu17 = +uy +uz;
const float cu18 = -uy -uz;
const float cu19 = -ux +uy -uz;
const float cu20 = +ux -uy +uz;
const float cu21 = -ux -uy +uz;
const float cu22 = +ux +uy -uz;
const float cu23 = +ux -uy -uz;
const float cu24 = -ux +uy +uz;
const float cu25 = -ux -uy -uz;
const float cu26 = +ux +uy +uz;

constexpr float w0  = 8.f/27.f;
constexpr float w1  = 2.f/27.f;
constexpr float w2  = 1.f/54.f;
constexpr float w3 = 1.f/216.f;

const float feq0  = rho*w0 *(1.f + 3.f*cu0  + 4.5f*cu0 *cu0  - 1.5f*u2);
const float feq1  = rho*w1 *(1.f + 3.f*cu1  + 4.5f*cu1 *cu1  - 1.5f*u2);
const float	feq2  = rho*w1 *(1.f + 3.f*cu2  + 4.5f*cu2 *cu2  - 1.5f*u2);
const float	feq3  = rho*w1 *(1.f + 3.f*cu3  + 4.5f*cu3 *cu3  - 1.5f*u2);
const float	feq4  = rho*w1 *(1.f + 3.f*cu4  + 4.5f*cu4 *cu4  - 1.5f*u2);
const float	feq5  = rho*w1 *(1.f + 3.f*cu5  + 4.5f*cu5 *cu5  - 1.5f*u2);
const float	feq6  = rho*w1 *(1.f + 3.f*cu6  + 4.5f*cu6 *cu6  - 1.5f*u2);
const float	feq7  = rho*w2 *(1.f + 3.f*cu7  + 4.5f*cu7 *cu7  - 1.5f*u2);
const float	feq8  = rho*w2 *(1.f + 3.f*cu8  + 4.5f*cu8 *cu8  - 1.5f*u2);
const float	feq9  = rho*w2 *(1.f + 3.f*cu9  + 4.5f*cu9 *cu9  - 1.5f*u2);
const float	feq10 = rho*w2 *(1.f + 3.f*cu10 + 4.5f*cu10*cu10 - 1.5f*u2);
const float	feq11 = rho*w2 *(1.f + 3.f*cu11 + 4.5f*cu11*cu11 - 1.5f*u2);
const float	feq12 = rho*w2 *(1.f + 3.f*cu12 + 4.5f*cu12*cu12 - 1.5f*u2);
const float	feq13 = rho*w2 *(1.f + 3.f*cu13 + 4.5f*cu13*cu13 - 1.5f*u2);
const float	feq14 = rho*w2 *(1.f + 3.f*cu14 + 4.5f*cu14*cu14 - 1.5f*u2);
const float	feq15 = rho*w2 *(1.f + 3.f*cu15 + 4.5f*cu15*cu15 - 1.5f*u2);
const float	feq16 = rho*w2 *(1.f + 3.f*cu16 + 4.5f*cu16*cu16 - 1.5f*u2);
const float	feq17 = rho*w2 *(1.f + 3.f*cu17 + 4.5f*cu17*cu17 - 1.5f*u2);
const float	feq18 = rho*w2 *(1.f + 3.f*cu18 + 4.5f*cu18*cu18 - 1.5f*u2);
const float	feq19 = rho*w3 *(1.f + 3.f*cu19 + 4.5f*cu19*cu19 - 1.5f*u2);
const float	feq20 = rho*w3 *(1.f + 3.f*cu20 + 4.5f*cu20*cu20 - 1.5f*u2);
const float	feq21 = rho*w3 *(1.f + 3.f*cu21 + 4.5f*cu21*cu21 - 1.5f*u2);
const float	feq22 = rho*w3 *(1.f + 3.f*cu22 + 4.5f*cu22*cu22 - 1.5f*u2);
const float	feq23 = rho*w3 *(1.f + 3.f*cu23 + 4.5f*cu23*cu23 - 1.5f*u2);
const float	feq24 = rho*w3 *(1.f + 3.f*cu24 + 4.5f*cu24*cu24 - 1.5f*u2);
const float	feq25 = rho*w3 *(1.f + 3.f*cu25 + 4.5f*cu25*cu25 - 1.5f*u2);
const float	feq26 = rho*w3 *(1.f + 3.f*cu26 + 4.5f*cu26*cu26 - 1.5f*u2);	
