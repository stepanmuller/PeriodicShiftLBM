// Perform half-way bounceback BC on f0, f1 ...

// id: { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
// cx: { 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
// cy: { 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
// cz: { 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };

float fTemp;

fTemp = f1;  f1  = f2;  f2  = fTemp;
fTemp = f3;  f3  = f4;  f4  = fTemp;
fTemp = f5;  f5  = f6;  f6  = fTemp;

fTemp = f7;  f7  = f8;  f8  = fTemp;
fTemp = f9;  f9  = f10; f10 = fTemp;
fTemp = f11; f11 = f12; f12 = fTemp;
fTemp = f13; f13 = f14; f14 = fTemp;
fTemp = f15; f15 = f16; f16 = fTemp;
fTemp = f17; f17 = f18; f18 = fTemp;

fTemp = f19; f19 = f20; f20 = fTemp;
fTemp = f21; f21 = f22; f22 = fTemp;
fTemp = f23; f23 = f24; f24 = fTemp;
fTemp = f25; f25 = f26; f26 = fTemp;
