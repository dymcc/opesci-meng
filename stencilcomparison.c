//Tiled Stencil
for (int time = 1; time < time_size - 1; time += 1)
{
    int skew = 4*time;
    for (xx=0;xx<=21;xx++) {
        for (yy=0;yy<=21;yy++) {
            for (zz=0;zz<=21;zz++) {
                for (x=max(4,24*xx);x<=24*xx+23;x++) {
                    for (y=max(4,24*yy);y<=24*yy+23;y++) {
                        for (z=max(4,24*zz);z<=24*zz+23;z++) {
                            float tcse0 = 3.04F*damp[x][y][z];
                        }
                    }
                }
            }
        }
    }
}

//Skewed & Tiled Stencil
for (int time = 1; time < time_size - 1; time += 1)
{
    int skew = 4*time;
    for (xx=0;xx<=21;xx++) {
        for (yy=0;yy<=21;yy++) {
            for (zz=0;zz<=21;zz++) {
                for (x=max(skew+4,skew+24*xx);x<=skew+24*xx+23;x++) {
                    for (y=max(skew+4,skew+24*yy);y<=skew+24*yy+23;y++) {
                        for (z=max(skew+4,skew+24*zz);z<=skew+24*zz+23;z++) {
                            float tcse0 = 3.04F*damp[x-skew][y-skew][z-skew];
                        }
                    }
                }
            }
        }
    }
}

//Skewed & Tiled Stencil with Interchange
for (xx=0;xx<=21;xx++) {
    for (yy=0;yy<=21;yy++) {
        for (zz=0;zz<=21;zz++) {
            for (int time = 1; time < time_size - 1; time += 1)
            {
                int skew = 4*time;
                for (x=max(skew+4,skew+24*xx);x<=skew+24*xx+23;x++) {
                    for (y=max(skew+4,skew+24*yy);y<=skew+24*yy+23;y++) {
                        for (z=max(skew+4,skew+24*zz);z<=skew+24*zz+23;z++) {
                            float tcse0 = 3.04F*damp[x-skew][y-skew][z-skew];
                        }
                    }
                }
            }
        }
    }
}

