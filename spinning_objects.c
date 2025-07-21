// spinning_objects.c
// ASCII 3D Shape Renderer in C
// Shapes: Donut, Sphere, Cube
// Press 'S' to stop the animation
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#ifdef _WIN32
#include <conio.h>
#include <windows.h>
#else
#include <unistd.h>
#include <termios.h>
#include <fcntl.h>
#endif

#define SCREEN_WIDTH 80
#define SCREEN_HEIGHT 24
#define THETA_SPACING 0.07
#define PHI_SPACING 0.02
#define R1 1
#define R2 2
#define K2 5
#define K1 (SCREEN_WIDTH * K2 * 3 / (8 * (R1 + R2)))

// Platform-specific sleep
void sleep_ms(int ms) {
#ifdef _WIN32
    Sleep(ms);
#else
    usleep(ms * 1000);
#endif
}

// Platform-specific kbhit and getch
#ifdef _WIN32
// Already included
#else
int kbhit(void) {
    struct termios oldt, newt;
    int ch;
    int oldf;
    tcgetattr(STDIN_FILENO, &oldt);
    newt = oldt;
    newt.c_lflag &= ~(ICANON | ECHO);
    tcsetattr(STDIN_FILENO, TCSANOW, &newt);
    oldf = fcntl(STDIN_FILENO, F_GETFL, 0);
    fcntl(STDIN_FILENO, F_SETFL, oldf | O_NONBLOCK);
    ch = getchar();
    tcsetattr(STDIN_FILENO, TCSANOW, &oldt);
    fcntl(STDIN_FILENO, F_SETFL, oldf);
    if(ch != EOF) {
        ungetc(ch, stdin);
        return 1;
    }
    return 0;
}
int getch(void) {
    struct termios oldt, newt;
    int ch;
    tcgetattr(STDIN_FILENO, &oldt);
    newt = oldt;
    newt.c_lflag &= ~(ICANON | ECHO);
    tcsetattr(STDIN_FILENO, TCSANOW, &newt);
    ch = getchar();
    tcsetattr(STDIN_FILENO, TCSANOW, &oldt);
    return ch;
}
#endif

// Donut animation (already high quality)
void animate_donut() {
    float A = 0, B = 0;
    char output[SCREEN_WIDTH * SCREEN_HEIGHT];
    float zbuffer[SCREEN_WIDTH * SCREEN_HEIGHT];
    const char *lum = ".,-~:;=!*#$@";
    int stop = 0;
    printf("\033[2J"); // Clear screen
    while (!stop) {
        memset(output, ' ', SCREEN_WIDTH * SCREEN_HEIGHT);
        memset(zbuffer, 0, sizeof(zbuffer));
        for (float theta = 0; theta < 2 * M_PI; theta += THETA_SPACING) {
            for (float phi = 0; phi < 2 * M_PI; phi += PHI_SPACING) {
                float sinA = sin(A), cosA = cos(A);
                float sinB = sin(B), cosB = cos(B);
                float sintheta = sin(theta), costheta = cos(theta);
                float sinphi = sin(phi), cosphi = cos(phi);
                float circlex = R2 + R1 * costheta;
                float circley = R1 * sintheta;
                float x = circlex * (cosB * cosphi + sinA * sinB * sinphi) - circley * cosA * sinB;
                float y = circlex * (sinB * cosphi - sinA * cosB * sinphi) + circley * cosA * cosB;
                float z = K2 + cosA * circlex * sinphi + circley * sinA;
                float ooz = 1 / z;
                int xp = (int)(SCREEN_WIDTH / 2 + K1 * ooz * x);
                int yp = (int)(SCREEN_HEIGHT / 2 - K1 * ooz * y);
                int idx = xp + yp * SCREEN_WIDTH;
                float L = cosphi * costheta * sinB - cosA * costheta * sinphi - sinA * sintheta + cosB * (cosA * sintheta - costheta * sinA * sinphi);
                if (idx >= 0 && idx < SCREEN_WIDTH * SCREEN_HEIGHT && L > 0) {
                    if (ooz > zbuffer[idx]) {
                        zbuffer[idx] = ooz;
                        int lum_idx = (int)(L * 8);
                        output[idx] = lum[lum_idx > 0 ? lum_idx : 0];
                    }
                }
            }
        }
        printf("\033[H");
        for (int i = 0; i < SCREEN_HEIGHT; ++i) {
            for (int j = 0; j < SCREEN_WIDTH; ++j) {
                putchar(output[i * SCREEN_WIDTH + j]);
            }
            putchar('\n');
        }
        A += 0.04;
        B += 0.02;
        sleep_ms(30);
#ifdef _WIN32
        if (_kbhit()) {
            char c = _getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#else
        if (kbhit()) {
            char c = getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#endif
    }
    printf("\nDonut animation stopped.\n");
}

// High-quality shaded sphere animation
void animate_sphere() {
    float A = 0, B = 0;
    char output[SCREEN_WIDTH * SCREEN_HEIGHT];
    float zbuffer[SCREEN_WIDTH * SCREEN_HEIGHT];
    const char *lum = ".,-~:;=!*#$@";
    int stop = 0;
    printf("\033[2J");
    while (!stop) {
        memset(output, ' ', SCREEN_WIDTH * SCREEN_HEIGHT);
        memset(zbuffer, 0, sizeof(zbuffer));
        for (float theta = 0; theta < 2 * M_PI; theta += THETA_SPACING) {
            for (float phi = 0; phi < 2 * M_PI; phi += PHI_SPACING) {
                // Sphere coordinates (radius 1)
                float sintheta = sin(theta), costheta = cos(theta);
                float sinphi = sin(phi), cosphi = cos(phi);
                float x = costheta * cosphi;
                float y = costheta * sinphi;
                float z = sintheta;
                // Rotate in 3D
                float x1 = x * cos(B) + z * sin(B);
                float y1 = y;
                float z1 = -x * sin(B) + z * cos(B);
                float x2 = x1;
                float y2 = y1 * cos(A) - z1 * sin(A);
                float z2 = y1 * sin(A) + z1 * cos(A);
                float K = 15;
                float ooz = 1 / (K2 + z2);
                int xp = (int)(SCREEN_WIDTH / 2 + K1 * ooz * x2);
                int yp = (int)(SCREEN_HEIGHT / 2 - K1 * ooz * y2);
                int idx = xp + yp * SCREEN_WIDTH;
                // Lighting: dot product of normal and light direction
                float Lx = 0, Ly = 1, Lz = -1; // Light direction
                float norm = sqrt(x * x + y * y + z * z);
                float L = (x * Lx + y * Ly + z * Lz) / (norm * sqrt(Lx * Lx + Ly * Ly + Lz * Lz));
                if (idx >= 0 && idx < SCREEN_WIDTH * SCREEN_HEIGHT && L > 0) {
                    if (ooz > zbuffer[idx]) {
                        zbuffer[idx] = ooz;
                        int lum_idx = (int)(L * 10);
                        if (lum_idx < 0) lum_idx = 0;
                        if (lum_idx > 11) lum_idx = 11;
                        output[idx] = lum[lum_idx];
                    }
                }
            }
        }
        printf("\033[H");
        for (int i = 0; i < SCREEN_HEIGHT; ++i) {
            for (int j = 0; j < SCREEN_WIDTH; ++j) {
                putchar(output[i * SCREEN_WIDTH + j]);
            }
            putchar('\n');
        }
        A += 0.04;
        B += 0.02;
        sleep_ms(30);
#ifdef _WIN32
        if (_kbhit()) {
            char c = _getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#else
        if (kbhit()) {
            char c = getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#endif
    }
    printf("\nSphere animation stopped.\n");
}

// High-quality shaded cube animation
void animate_cube() {
    float A = 0, B = 0;
    char output[SCREEN_WIDTH * SCREEN_HEIGHT];
    float zbuffer[SCREEN_WIDTH * SCREEN_HEIGHT];
    const char *lum = ".,-~:;=!*#$@";
    int stop = 0;
    printf("\033[2J");
    // For each face, sample points densely
    // Cube faces: 6 faces, each face is a square in 3D
    float faces[6][4][3] = {
        // +X
        {{1,1,1},{1,1,-1},{1,-1,-1},{1,-1,1}},
        // -X
        {{-1,1,1},{-1,1,-1},{-1,-1,-1},{-1,-1,1}},
        // +Y
        {{-1,1,1},{1,1,1},{1,1,-1},{-1,1,-1}},
        // -Y
        {{-1,-1,1},{1,-1,1},{1,-1,-1},{-1,-1,-1}},
        // +Z
        {{-1,1,1},{1,1,1},{1,-1,1},{-1,-1,1}},
        // -Z
        {{-1,1,-1},{1,1,-1},{1,-1,-1},{-1,-1,-1}}
    };
    float face_normals[6][3] = {
        {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1}
    };
    while (!stop) {
        memset(output, ' ', SCREEN_WIDTH * SCREEN_HEIGHT);
        memset(zbuffer, 0, sizeof(zbuffer));
        for (int f = 0; f < 6; ++f) {
            // Sample each face densely
            for (float u = 0; u <= 1.0; u += 0.03) {
                for (float v = 0; v <= 1.0; v += 0.03) {
                    // Bilinear interpolation for face
                    float p[3] = {
                        (1-u)*(1-v)*faces[f][0][0] + u*(1-v)*faces[f][1][0] + u*v*faces[f][2][0] + (1-u)*v*faces[f][3][0],
                        (1-u)*(1-v)*faces[f][0][1] + u*(1-v)*faces[f][1][1] + u*v*faces[f][2][1] + (1-u)*v*faces[f][3][1],
                        (1-u)*(1-v)*faces[f][0][2] + u*(1-v)*faces[f][1][2] + u*v*faces[f][2][2] + (1-u)*v*faces[f][3][2]
                    };
                    // Rotate in 3D
                    float x = p[0], y = p[1], z = p[2];
                    float x1 = x * cos(B) + z * sin(B);
                    float y1 = y;
                    float z1 = -x * sin(B) + z * cos(B);
                    float x2 = x1;
                    float y2 = y1 * cos(A) - z1 * sin(A);
                    float z2 = y1 * sin(A) + z1 * cos(A);
                    float K = 8;
                    float ooz = 1 / (K2 + z2);
                    int xp = (int)(SCREEN_WIDTH / 2 + K1 * ooz * x2);
                    int yp = (int)(SCREEN_HEIGHT / 2 - K1 * ooz * y2);
                    int idx = xp + yp * SCREEN_WIDTH;
                    // Lighting: dot product of face normal and light direction
                    float nx = face_normals[f][0], ny = face_normals[f][1], nz = face_normals[f][2];
                    // Rotate normal
                    float nx1 = nx * cos(B) + nz * sin(B);
                    float ny1 = ny;
                    float nz1 = -nx * sin(B) + nz * cos(B);
                    float nx2 = nx1;
                    float ny2 = ny1 * cos(A) - nz1 * sin(A);
                    float nz2 = ny1 * sin(A) + nz1 * cos(A);
                    float Lx = 0, Ly = 1, Lz = -1; // Light direction
                    float norm = sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);
                    float L = (nx2 * Lx + ny2 * Ly + nz2 * Lz) / (norm * sqrt(Lx * Lx + Ly * Ly + Lz * Lz));
                    if (idx >= 0 && idx < SCREEN_WIDTH * SCREEN_HEIGHT && L > 0) {
                        if (ooz > zbuffer[idx]) {
                            zbuffer[idx] = ooz;
                            int lum_idx = (int)(L * 10);
                            if (lum_idx < 0) lum_idx = 0;
                            if (lum_idx > 11) lum_idx = 11;
                            output[idx] = lum[lum_idx];
                        }
                    }
                }
            }
        }
        printf("\033[H");
        for (int i = 0; i < SCREEN_HEIGHT; ++i) {
            for (int j = 0; j < SCREEN_WIDTH; ++j) {
                putchar(output[i * SCREEN_WIDTH + j]);
            }
            putchar('\n');
        }
        A += 0.04;
        B += 0.02;
        sleep_ms(30);
#ifdef _WIN32
        if (_kbhit()) {
            char c = _getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#else
        if (kbhit()) {
            char c = getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#endif
    }
    printf("\nCube animation stopped.\n");
}

// High-quality shaded pyramid animation
void animate_pyramid() {
    float A = 0, B = 0;
    char output[SCREEN_WIDTH * SCREEN_HEIGHT];
    float zbuffer[SCREEN_WIDTH * SCREEN_HEIGHT];
    const char *lum = ".,-~:;=!*#$@";
    int stop = 0;
    printf("\033[2J");
    // Pyramid: square base, 4 triangular faces, centered at origin, height 2
    float base[4][3] = {
        {-1, -1, -1}, {1, -1, -1}, {1, -1, 1}, {-1, -1, 1}
    };
    float apex[3] = {0, 1, 0};
    // Each face: 3 vertices (apex, base[i], base[(i+1)%4])
    while (!stop) {
        memset(output, ' ', SCREEN_WIDTH * SCREEN_HEIGHT);
        memset(zbuffer, 0, sizeof(zbuffer));
        // Draw base (as a filled quad)
        for (float u = 0; u <= 1.0; u += 0.03) {
            for (float v = 0; v <= 1.0 - u; v += 0.03) {
                float w = 1.0 - u - v;
                float x = w * base[0][0] + u * base[1][0] + v * base[2][0];
                float y = w * base[0][1] + u * base[1][1] + v * base[2][1];
                float z = w * base[0][2] + u * base[1][2] + v * base[2][2];
                // Normal for base: (0, -1, 0)
                float nx = 0, ny = -1, nz = 0;
                // Rotate point and normal
                float x1 = x * cos(B) + z * sin(B);
                float y1 = y;
                float z1 = -x * sin(B) + z * cos(B);
                float x2 = x1;
                float y2 = y1 * cos(A) - z1 * sin(A);
                float z2 = y1 * sin(A) + z1 * cos(A);
                float nx1 = nx * cos(B) + nz * sin(B);
                float ny1 = ny;
                float nz1 = -nx * sin(B) + nz * cos(B);
                float nx2 = nx1;
                float ny2 = ny1 * cos(A) - nz1 * sin(A);
                float nz2 = ny1 * sin(A) + nz1 * cos(A);
                float K = 8;
                float ooz = 1 / (K2 + z2);
                int xp = (int)(SCREEN_WIDTH / 2 + K1 * ooz * x2);
                int yp = (int)(SCREEN_HEIGHT / 2 - K1 * ooz * y2);
                int idx = xp + yp * SCREEN_WIDTH;
                float Lx = 0, Ly = 1, Lz = -1;
                float norm = sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);
                float L = (nx2 * Lx + ny2 * Ly + nz2 * Lz) / (norm * sqrt(Lx * Lx + Ly * Ly + Lz * Lz));
                if (idx >= 0 && idx < SCREEN_WIDTH * SCREEN_HEIGHT && L > 0) {
                    if (ooz > zbuffer[idx]) {
                        zbuffer[idx] = ooz;
                        int lum_idx = (int)(L * 10);
                        if (lum_idx < 0) lum_idx = 0;
                        if (lum_idx > 11) lum_idx = 11;
                        output[idx] = lum[lum_idx];
                    }
                }
            }
        }
        // Draw 4 triangular faces
        for (int f = 0; f < 4; ++f) {
            float *v0 = apex;
            float *v1 = base[f];
            float *v2 = base[(f+1)%4];
            for (float u = 0; u <= 1.0; u += 0.03) {
                for (float v = 0; v <= 1.0 - u; v += 0.03) {
                    float w = 1.0 - u - v;
                    float x = w * v0[0] + u * v1[0] + v * v2[0];
                    float y = w * v0[1] + u * v1[1] + v * v2[1];
                    float z = w * v0[2] + u * v1[2] + v * v2[2];
                    // Normal for face
                    float ux = v1[0] - v0[0], uy = v1[1] - v0[1], uz = v1[2] - v0[2];
                    float vx_ = v2[0] - v0[0], vy_ = v2[1] - v0[1], vz_ = v2[2] - v0[2];
                    float nx = uy * vz_ - uz * vy_;
                    float ny = uz * vx_ - ux * vz_;
                    float nz = ux * vy_ - uy * vx_;
                    // Rotate point and normal
                    float x1 = x * cos(B) + z * sin(B);
                    float y1 = y;
                    float z1 = -x * sin(B) + z * cos(B);
                    float x2 = x1;
                    float y2 = y1 * cos(A) - z1 * sin(A);
                    float z2 = y1 * sin(A) + z1 * cos(A);
                    float nx1 = nx * cos(B) + nz * sin(B);
                    float ny1 = ny;
                    float nz1 = -nx * sin(B) + nz * cos(B);
                    float nx2 = nx1;
                    float ny2 = ny1 * cos(A) - nz1 * sin(A);
                    float nz2 = ny1 * sin(A) + nz1 * cos(A);
                    float K = 8;
                    float ooz = 1 / (K2 + z2);
                    int xp = (int)(SCREEN_WIDTH / 2 + K1 * ooz * x2);
                    int yp = (int)(SCREEN_HEIGHT / 2 - K1 * ooz * y2);
                    int idx = xp + yp * SCREEN_WIDTH;
                    float Lx = 0, Ly = 1, Lz = -1;
                    float norm = sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);
                    float L = (nx2 * Lx + ny2 * Ly + nz2 * Lz) / (norm * sqrt(Lx * Lx + Ly * Ly + Lz * Lz));
                    if (idx >= 0 && idx < SCREEN_WIDTH * SCREEN_HEIGHT && L > 0) {
                        if (ooz > zbuffer[idx]) {
                            zbuffer[idx] = ooz;
                            int lum_idx = (int)(L * 10);
                            if (lum_idx < 0) lum_idx = 0;
                            if (lum_idx > 11) lum_idx = 11;
                            output[idx] = lum[lum_idx];
                        }
                    }
                }
            }
        }
        printf("\033[H");
        for (int i = 0; i < SCREEN_HEIGHT; ++i) {
            for (int j = 0; j < SCREEN_WIDTH; ++j) {
                putchar(output[i * SCREEN_WIDTH + j]);
            }
            putchar('\n');
        }
        A += 0.04;
        B += 0.02;
        sleep_ms(30);
#ifdef _WIN32
        if (_kbhit()) {
            char c = _getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#else
        if (kbhit()) {
            char c = getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#endif
    }
    printf("\nPyramid animation stopped.\n");
}

// High-quality shaded cone animation
void animate_cone() {
    float A = 0, B = 0;
    char output[SCREEN_WIDTH * SCREEN_HEIGHT];
    float zbuffer[SCREEN_WIDTH * SCREEN_HEIGHT];
    const char *lum = ".,-~:;=!*#$@";
    int stop = 0;
    printf("\033[2J");
    // Cone: base at y = -1, apex at y = 1, radius 1
    while (!stop) {
        memset(output, ' ', SCREEN_WIDTH * SCREEN_HEIGHT);
        memset(zbuffer, 0, sizeof(zbuffer));
        // Draw base (disk)
        for (float r = 0; r <= 1.0; r += 0.03) {
            for (float phi = 0; phi < 2 * M_PI; phi += PHI_SPACING) {
                float x = r * cos(phi);
                float y = -1;
                float z = r * sin(phi);
                // Normal for base: (0, -1, 0)
                float nx = 0, ny = -1, nz = 0;
                // Rotate point and normal
                float x1 = x * cos(B) + z * sin(B);
                float y1 = y;
                float z1 = -x * sin(B) + z * cos(B);
                float x2 = x1;
                float y2 = y1 * cos(A) - z1 * sin(A);
                float z2 = y1 * sin(A) + z1 * cos(A);
                float nx1 = nx * cos(B) + nz * sin(B);
                float ny1 = ny;
                float nz1 = -nx * sin(B) + nz * cos(B);
                float nx2 = nx1;
                float ny2 = ny1 * cos(A) - nz1 * sin(A);
                float nz2 = ny1 * sin(A) + nz1 * cos(A);
                float K = 8;
                float ooz = 1 / (K2 + z2);
                int xp = (int)(SCREEN_WIDTH / 2 + K1 * ooz * x2);
                int yp = (int)(SCREEN_HEIGHT / 2 - K1 * ooz * y2);
                int idx = xp + yp * SCREEN_WIDTH;
                float Lx = 0, Ly = 1, Lz = -1;
                float norm = sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);
                float L = (nx2 * Lx + ny2 * Ly + nz2 * Lz) / (norm * sqrt(Lx * Lx + Ly * Ly + Lz * Lz));
                if (idx >= 0 && idx < SCREEN_WIDTH * SCREEN_HEIGHT && L > 0) {
                    if (ooz > zbuffer[idx]) {
                        zbuffer[idx] = ooz;
                        int lum_idx = (int)(L * 10);
                        if (lum_idx < 0) lum_idx = 0;
                        if (lum_idx > 11) lum_idx = 11;
                        output[idx] = lum[lum_idx];
                    }
                }
            }
        }
        // Draw cone surface
        for (float h = 0; h <= 1.0; h += 0.02) {
            float y = -1 + 2 * h;
            float r = 1 - h;
            for (float phi = 0; phi < 2 * M_PI; phi += PHI_SPACING) {
                float x = r * cos(phi);
                float z = r * sin(phi);
                // Apex at (0,1,0), base at y=-1
                float v0[3] = {0, 1, 0};
                float v1[3] = {x, y, z};
                float v2[3] = {r * cos(phi + PHI_SPACING), y, r * sin(phi + PHI_SPACING)};
                // Normal for surface triangle
                float ux = v1[0] - v0[0], uy = v1[1] - v0[1], uz = v1[2] - v0[2];
                float vx_ = v2[0] - v0[0], vy_ = v2[1] - v0[1], vz_ = v2[2] - v0[2];
                float nx = uy * vz_ - uz * vy_;
                float ny = uz * vx_ - ux * vz_;
                float nz = ux * vy_ - uy * vx_;
                // Point on surface (midpoint of triangle)
                float px = (v0[0] + v1[0] + v2[0]) / 3.0;
                float py = (v0[1] + v1[1] + v2[1]) / 3.0;
                float pz = (v0[2] + v1[2] + v2[2]) / 3.0;
                // Rotate point and normal
                float x1 = px * cos(B) + pz * sin(B);
                float y1 = py;
                float z1 = -px * sin(B) + pz * cos(B);
                float x2 = x1;
                float y2 = y1 * cos(A) - z1 * sin(A);
                float z2 = y1 * sin(A) + z1 * cos(A);
                float nx1 = nx * cos(B) + nz * sin(B);
                float ny1 = ny;
                float nz1 = -nx * sin(B) + nz * cos(B);
                float nx2 = nx1;
                float ny2 = ny1 * cos(A) - nz1 * sin(A);
                float nz2 = ny1 * sin(A) + nz1 * cos(A);
                float K = 8;
                float ooz = 1 / (K2 + z2);
                int xp = (int)(SCREEN_WIDTH / 2 + K1 * ooz * x2);
                int yp = (int)(SCREEN_HEIGHT / 2 - K1 * ooz * y2);
                int idx = xp + yp * SCREEN_WIDTH;
                float Lx = 0, Ly = 1, Lz = -1;
                float norm = sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);
                float L = (nx2 * Lx + ny2 * Ly + nz2 * Lz) / (norm * sqrt(Lx * Lx + Ly * Ly + Lz * Lz));
                if (idx >= 0 && idx < SCREEN_WIDTH * SCREEN_HEIGHT && L > 0) {
                    if (ooz > zbuffer[idx]) {
                        zbuffer[idx] = ooz;
                        int lum_idx = (int)(L * 10);
                        if (lum_idx < 0) lum_idx = 0;
                        if (lum_idx > 11) lum_idx = 11;
                        output[idx] = lum[lum_idx];
                    }
                }
            }
        }
        printf("\033[H");
        for (int i = 0; i < SCREEN_HEIGHT; ++i) {
            for (int j = 0; j < SCREEN_WIDTH; ++j) {
                putchar(output[i * SCREEN_WIDTH + j]);
            }
            putchar('\n');
        }
        A += 0.04;
        B += 0.02;
        sleep_ms(30);
#ifdef _WIN32
        if (_kbhit()) {
            char c = _getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#else
        if (kbhit()) {
            char c = getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#endif
    }
    printf("\nCone animation stopped.\n");
}

// High-quality shaded cylinder animation
void animate_cylinder() {
    float A = 0, B = 0;
    char output[SCREEN_WIDTH * SCREEN_HEIGHT];
    float zbuffer[SCREEN_WIDTH * SCREEN_HEIGHT];
    const char *lum = ".,-~:;=!*#$@";
    int stop = 0;
    printf("\033[2J");
    // Cylinder: height 2 (y from -1 to 1), radius 1
    while (!stop) {
        memset(output, ' ', SCREEN_WIDTH * SCREEN_HEIGHT);
        memset(zbuffer, 0, sizeof(zbuffer));
        // Draw top and bottom disks
        for (int cap = -1; cap <= 1; cap += 2) {
            for (float r = 0; r <= 1.0; r += 0.03) {
                for (float phi = 0; phi < 2 * M_PI; phi += PHI_SPACING) {
                    float x = r * cos(phi);
                    float y = cap;
                    float z = r * sin(phi);
                    // Normal for cap: (0, cap, 0)
                    float nx = 0, ny = cap, nz = 0;
                    // Rotate point and normal
                    float x1 = x * cos(B) + z * sin(B);
                    float y1 = y;
                    float z1 = -x * sin(B) + z * cos(B);
                    float x2 = x1;
                    float y2 = y1 * cos(A) - z1 * sin(A);
                    float z2 = y1 * sin(A) + z1 * cos(A);
                    float nx1 = nx * cos(B) + nz * sin(B);
                    float ny1 = ny;
                    float nz1 = -nx * sin(B) + nz * cos(B);
                    float nx2 = nx1;
                    float ny2 = ny1 * cos(A) - nz1 * sin(A);
                    float nz2 = ny1 * sin(A) + nz1 * cos(A);
                    float K = 8;
                    float ooz = 1 / (K2 + z2);
                    int xp = (int)(SCREEN_WIDTH / 2 + K1 * ooz * x2);
                    int yp = (int)(SCREEN_HEIGHT / 2 - K1 * ooz * y2);
                    int idx = xp + yp * SCREEN_WIDTH;
                    float Lx = 0, Ly = 1, Lz = -1;
                    float norm = sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);
                    float L = (nx2 * Lx + ny2 * Ly + nz2 * Lz) / (norm * sqrt(Lx * Lx + Ly * Ly + Lz * Lz));
                    if (idx >= 0 && idx < SCREEN_WIDTH * SCREEN_HEIGHT && L > 0) {
                        if (ooz > zbuffer[idx]) {
                            zbuffer[idx] = ooz;
                            int lum_idx = (int)(L * 10);
                            if (lum_idx < 0) lum_idx = 0;
                            if (lum_idx > 11) lum_idx = 11;
                            output[idx] = lum[lum_idx];
                        }
                    }
                }
            }
        }
        // Draw side surface
        for (float y = -1; y <= 1.0; y += 0.03) {
            for (float phi = 0; phi < 2 * M_PI; phi += PHI_SPACING) {
                float x = cos(phi);
                float z = sin(phi);
                // Normal for side: (cos(phi), 0, sin(phi))
                float nx = cos(phi), ny = 0, nz = sin(phi);
                // Rotate point and normal
                float x1 = x * cos(B) + z * sin(B);
                float y1 = y;
                float z1 = -x * sin(B) + z * cos(B);
                float x2 = x1;
                float y2 = y1 * cos(A) - z1 * sin(A);
                float z2 = y1 * sin(A) + z1 * cos(A);
                float nx1 = nx * cos(B) + nz * sin(B);
                float ny1 = ny;
                float nz1 = -nx * sin(B) + nz * cos(B);
                float nx2 = nx1;
                float ny2 = ny1 * cos(A) - nz1 * sin(A);
                float nz2 = ny1 * sin(A) + nz1 * cos(A);
                float K = 8;
                float ooz = 1 / (K2 + z2);
                int xp = (int)(SCREEN_WIDTH / 2 + K1 * ooz * x2);
                int yp = (int)(SCREEN_HEIGHT / 2 - K1 * ooz * y2);
                int idx = xp + yp * SCREEN_WIDTH;
                float Lx = 0, Ly = 1, Lz = -1;
                float norm = sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);
                float L = (nx2 * Lx + ny2 * Ly + nz2 * Lz) / (norm * sqrt(Lx * Lx + Ly * Ly + Lz * Lz));
                if (idx >= 0 && idx < SCREEN_WIDTH * SCREEN_HEIGHT && L > 0) {
                    if (ooz > zbuffer[idx]) {
                        zbuffer[idx] = ooz;
                        int lum_idx = (int)(L * 10);
                        if (lum_idx < 0) lum_idx = 0;
                        if (lum_idx > 11) lum_idx = 11;
                        output[idx] = lum[lum_idx];
                    }
                }
            }
        }
        printf("\033[H");
        for (int i = 0; i < SCREEN_HEIGHT; ++i) {
            for (int j = 0; j < SCREEN_WIDTH; ++j) {
                putchar(output[i * SCREEN_WIDTH + j]);
            }
            putchar('\n');
        }
        A += 0.04;
        B += 0.02;
        sleep_ms(30);
#ifdef _WIN32
        if (_kbhit()) {
            char c = _getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#else
        if (kbhit()) {
            char c = getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#endif
    }
    printf("\nCylinder animation stopped.\n");
}

// High-quality shaded hourglass animation
void animate_hourglass() {
    float A = 0, B = 0;
    char output[SCREEN_WIDTH * SCREEN_HEIGHT];
    float zbuffer[SCREEN_WIDTH * SCREEN_HEIGHT];
    const char *lum = ".,-~:;=!*#$@";
    int stop = 0;
    printf("\033[2J");
    // Hourglass: two cones joined at y=0, base radius 1, height 2 (y from -1 to 1)
    while (!stop) {
        memset(output, ' ', SCREEN_WIDTH * SCREEN_HEIGHT);
        memset(zbuffer, 0, sizeof(zbuffer));
        // Top cone (apex at y=1, base at y=0)
        for (float h = 0; h <= 1.0; h += 0.02) {
            float y = 1 - h;
            float r = h;
            for (float phi = 0; phi < 2 * M_PI; phi += PHI_SPACING) {
                float x = r * cos(phi);
                float z = r * sin(phi);
                // Apex at (0,1,0), base at y=0
                float v0[3] = {0, 1, 0};
                float v1[3] = {x, y, z};
                float v2[3] = {r * cos(phi + PHI_SPACING), y, r * sin(phi + PHI_SPACING)};
                // Normal for surface triangle
                float ux = v1[0] - v0[0], uy = v1[1] - v0[1], uz = v1[2] - v0[2];
                float vx_ = v2[0] - v0[0], vy_ = v2[1] - v0[1], vz_ = v2[2] - v0[2];
                float nx = uy * vz_ - uz * vy_;
                float ny = uz * vx_ - ux * vz_;
                float nz = ux * vy_ - uy * vx_;
                // Point on surface (midpoint of triangle)
                float px = (v0[0] + v1[0] + v2[0]) / 3.0;
                float py = (v0[1] + v1[1] + v2[1]) / 3.0;
                float pz = (v0[2] + v1[2] + v2[2]) / 3.0;
                // Rotate point and normal
                float x1 = px * cos(B) + pz * sin(B);
                float y1 = py;
                float z1 = -px * sin(B) + pz * cos(B);
                float x2 = x1;
                float y2 = y1 * cos(A) - z1 * sin(A);
                float z2 = y1 * sin(A) + z1 * cos(A);
                float nx1 = nx * cos(B) + nz * sin(B);
                float ny1 = ny;
                float nz1 = -nx * sin(B) + nz * cos(B);
                float nx2 = nx1;
                float ny2 = ny1 * cos(A) - nz1 * sin(A);
                float nz2 = ny1 * sin(A) + nz1 * cos(A);
                float K = 8;
                float ooz = 1 / (K2 + z2);
                int xp = (int)(SCREEN_WIDTH / 2 + K1 * ooz * x2);
                int yp = (int)(SCREEN_HEIGHT / 2 - K1 * ooz * y2);
                int idx = xp + yp * SCREEN_WIDTH;
                float Lx = 0, Ly = 1, Lz = -1;
                float norm = sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);
                float L = (nx2 * Lx + ny2 * Ly + nz2 * Lz) / (norm * sqrt(Lx * Lx + Ly * Ly + Lz * Lz));
                if (idx >= 0 && idx < SCREEN_WIDTH * SCREEN_HEIGHT && L > 0) {
                    if (ooz > zbuffer[idx]) {
                        zbuffer[idx] = ooz;
                        int lum_idx = (int)(L * 10);
                        if (lum_idx < 0) lum_idx = 0;
                        if (lum_idx > 11) lum_idx = 11;
                        output[idx] = lum[lum_idx];
                    }
                }
            }
        }
        // Bottom cone (apex at y=-1, base at y=0)
        for (float h = 0; h <= 1.0; h += 0.02) {
            float y = -1 + h;
            float r = h;
            for (float phi = 0; phi < 2 * M_PI; phi += PHI_SPACING) {
                float x = r * cos(phi);
                float z = r * sin(phi);
                // Apex at (0,-1,0), base at y=0
                float v0[3] = {0, -1, 0};
                float v1[3] = {x, y, z};
                float v2[3] = {r * cos(phi + PHI_SPACING), y, r * sin(phi + PHI_SPACING)};
                // Normal for surface triangle
                float ux = v1[0] - v0[0], uy = v1[1] - v0[1], uz = v1[2] - v0[2];
                float vx_ = v2[0] - v0[0], vy_ = v2[1] - v0[1], vz_ = v2[2] - v0[2];
                float nx = uy * vz_ - uz * vy_;
                float ny = uz * vx_ - ux * vz_;
                float nz = ux * vy_ - uy * vx_;
                // Point on surface (midpoint of triangle)
                float px = (v0[0] + v1[0] + v2[0]) / 3.0;
                float py = (v0[1] + v1[1] + v2[1]) / 3.0;
                float pz = (v0[2] + v1[2] + v2[2]) / 3.0;
                // Rotate point and normal
                float x1 = px * cos(B) + pz * sin(B);
                float y1 = py;
                float z1 = -px * sin(B) + pz * cos(B);
                float x2 = x1;
                float y2 = y1 * cos(A) - z1 * sin(A);
                float z2 = y1 * sin(A) + z1 * cos(A);
                float nx1 = nx * cos(B) + nz * sin(B);
                float ny1 = ny;
                float nz1 = -nx * sin(B) + nz * cos(B);
                float nx2 = nx1;
                float ny2 = ny1 * cos(A) - nz1 * sin(A);
                float nz2 = ny1 * sin(A) + nz1 * cos(A);
                float K = 8;
                float ooz = 1 / (K2 + z2);
                int xp = (int)(SCREEN_WIDTH / 2 + K1 * ooz * x2);
                int yp = (int)(SCREEN_HEIGHT / 2 - K1 * ooz * y2);
                int idx = xp + yp * SCREEN_WIDTH;
                float Lx = 0, Ly = 1, Lz = -1;
                float norm = sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);
                float L = (nx2 * Lx + ny2 * Ly + nz2 * Lz) / (norm * sqrt(Lx * Lx + Ly * Ly + Lz * Lz));
                if (idx >= 0 && idx < SCREEN_WIDTH * SCREEN_HEIGHT && L > 0) {
                    if (ooz > zbuffer[idx]) {
                        zbuffer[idx] = ooz;
                        int lum_idx = (int)(L * 10);
                        if (lum_idx < 0) lum_idx = 0;
                        if (lum_idx > 11) lum_idx = 11;
                        output[idx] = lum[lum_idx];
                    }
                }
            }
        }
        printf("\033[H");
        for (int i = 0; i < SCREEN_HEIGHT; ++i) {
            for (int j = 0; j < SCREEN_WIDTH; ++j) {
                putchar(output[i * SCREEN_WIDTH + j]);
            }
            putchar('\n');
        }
        A += 0.04;
        B += 0.02;
        sleep_ms(30);
#ifdef _WIN32
        if (_kbhit()) {
            char c = _getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#else
        if (kbhit()) {
            char c = getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#endif
    }
    printf("\nHourglass animation stopped.\n");
}

// High-quality shaded icosahedron animation
void animate_icosahedron() {
    float A = 0, B = 0;
    char output[SCREEN_WIDTH * SCREEN_HEIGHT];
    float zbuffer[SCREEN_WIDTH * SCREEN_HEIGHT];
    const char *lum = ".,-~:;=!*#$@";
    int stop = 0;
    printf("\033[2J");
    // Icosahedron vertices and faces
    const float X = 0.525731112119133606f;
    const float Z = 0.850650808352039932f;
    float v[12][3] = {
        {-X, 0, Z}, {X, 0, Z}, {-X, 0, -Z}, {X, 0, -Z},
        {0, Z, X}, {0, Z, -X}, {0, -Z, X}, {0, -Z, -X},
        {Z, X, 0}, {-Z, X, 0}, {Z, -X, 0}, {-Z, -X, 0}
    };
    int f[20][3] = {
        {0,4,1},{0,9,4},{9,5,4},{4,5,8},{4,8,1},
        {8,10,1},{8,3,10},{5,3,8},{5,2,3},{2,7,3},
        {7,10,3},{7,6,10},{7,11,6},{11,0,6},{0,1,6},
        {6,1,10},{9,0,11},{9,11,2},{9,2,5},{7,2,11}
    };
    while (!stop) {
        memset(output, ' ', SCREEN_WIDTH * SCREEN_HEIGHT);
        memset(zbuffer, 0, sizeof(zbuffer));
        for (int i = 0; i < 20; ++i) {
            float *v0 = v[f[i][0]];
            float *v1 = v[f[i][1]];
            float *v2 = v[f[i][2]];
            for (float u = 0; u <= 1.0; u += 0.04) {
                for (float v_ = 0; v_ <= 1.0 - u; v_ += 0.04) {
                    float w = 1.0 - u - v_;
                    float x = w * v0[0] + u * v1[0] + v_ * v2[0];
                    float y = w * v0[1] + u * v1[1] + v_ * v2[1];
                    float z = w * v0[2] + u * v1[2] + v_ * v2[2];
                    // Normal for face
                    float ux = v1[0] - v0[0], uy = v1[1] - v0[1], uz = v1[2] - v0[2];
                    float vx_ = v2[0] - v0[0], vy_ = v2[1] - v0[1], vz_ = v2[2] - v0[2];
                    float nx = uy * vz_ - uz * vy_;
                    float ny = uz * vx_ - ux * vz_;
                    float nz = ux * vy_ - uy * vx_;
                    // Rotate point and normal
                    float x1 = x * cos(B) + z * sin(B);
                    float y1 = y;
                    float z1 = -x * sin(B) + z * cos(B);
                    float x2 = x1;
                    float y2 = y1 * cos(A) - z1 * sin(A);
                    float z2 = y1 * sin(A) + z1 * cos(A);
                    float nx1 = nx * cos(B) + nz * sin(B);
                    float ny1 = ny;
                    float nz1 = -nx * sin(B) + nz * cos(B);
                    float nx2 = nx1;
                    float ny2 = ny1 * cos(A) - nz1 * sin(A);
                    float nz2 = ny1 * sin(A) + nz1 * cos(A);
                    float K = 8;
                    float ooz = 1 / (K2 + z2);
                    int xp = (int)(SCREEN_WIDTH / 2 + K1 * ooz * x2);
                    int yp = (int)(SCREEN_HEIGHT / 2 - K1 * ooz * y2);
                    int idx = xp + yp * SCREEN_WIDTH;
                    float Lx = 0, Ly = 1, Lz = -1;
                    float norm = sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);
                    float L = (nx2 * Lx + ny2 * Ly + nz2 * Lz) / (norm * sqrt(Lx * Lx + Ly * Ly + Lz * Lz));
                    if (idx >= 0 && idx < SCREEN_WIDTH * SCREEN_HEIGHT && L > 0) {
                        if (ooz > zbuffer[idx]) {
                            zbuffer[idx] = ooz;
                            int lum_idx = (int)(L * 10);
                            if (lum_idx < 0) lum_idx = 0;
                            if (lum_idx > 11) lum_idx = 11;
                            output[idx] = lum[lum_idx];
                        }
                    }
                }
            }
        }
        printf("\033[H");
        for (int i = 0; i < SCREEN_HEIGHT; ++i) {
            for (int j = 0; j < SCREEN_WIDTH; ++j) {
                putchar(output[i * SCREEN_WIDTH + j]);
            }
            putchar('\n');
        }
        A += 0.04;
        B += 0.02;
        sleep_ms(30);
#ifdef _WIN32
        if (_kbhit()) {
            char c = _getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#else
        if (kbhit()) {
            char c = getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#endif
    }
    printf("\nIcosahedron animation stopped.\n");
}

// High-quality shaded buckyball (truncated icosahedron) animation
// For simplicity, render as a set of small triangles on a sphere (approximation)
void animate_buckyball() {
    float A = 0, B = 0;
    char output[SCREEN_WIDTH * SCREEN_HEIGHT];
    float zbuffer[SCREEN_WIDTH * SCREEN_HEIGHT];
    const char *lum = ".,-~:;=!*#$@";
    int stop = 0;
    printf("\033[2J");
    // Approximate buckyball by subdividing a sphere into many small triangles
    float radius = 1.0;
    for (;!stop;) {
        memset(output, ' ', SCREEN_WIDTH * SCREEN_HEIGHT);
        memset(zbuffer, 0, sizeof(zbuffer));
        for (float theta = 0; theta < M_PI; theta += 0.15) {
            for (float phi = 0; phi < 2 * M_PI; phi += 0.15) {
                // Center of triangle patch
                float sintheta = sin(theta), costheta = cos(theta);
                float sinphi = sin(phi), cosphi = cos(phi);
                float x = radius * costheta * cosphi;
                float y = radius * costheta * sinphi;
                float z = radius * sintheta;
                // Normal is just the position vector (since it's a sphere)
                float nx = x, ny = y, nz = z;
                // Rotate point and normal
                float x1 = x * cos(B) + z * sin(B);
                float y1 = y;
                float z1 = -x * sin(B) + z * cos(B);
                float x2 = x1;
                float y2 = y1 * cos(A) - z1 * sin(A);
                float z2 = y1 * sin(A) + z1 * cos(A);
                float nx1 = nx * cos(B) + nz * sin(B);
                float ny1 = ny;
                float nz1 = -nx * sin(B) + nz * cos(B);
                float nx2 = nx1;
                float ny2 = ny1 * cos(A) - nz1 * sin(A);
                float nz2 = ny1 * sin(A) + nz1 * cos(A);
                float K = 15;
                float ooz = 1 / (K2 + z2);
                int xp = (int)(SCREEN_WIDTH / 2 + K1 * ooz * x2);
                int yp = (int)(SCREEN_HEIGHT / 2 - K1 * ooz * y2);
                int idx = xp + yp * SCREEN_WIDTH;
                float Lx = 0, Ly = 1, Lz = -1;
                float norm = sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);
                float L = (nx2 * Lx + ny2 * Ly + nz2 * Lz) / (norm * sqrt(Lx * Lx + Ly * Ly + Lz * Lz));
                if (idx >= 0 && idx < SCREEN_WIDTH * SCREEN_HEIGHT && L > 0) {
                    if (ooz > zbuffer[idx]) {
                        zbuffer[idx] = ooz;
                        int lum_idx = (int)(L * 10);
                        if (lum_idx < 0) lum_idx = 0;
                        if (lum_idx > 11) lum_idx = 11;
                        output[idx] = lum[lum_idx];
                    }
                }
            }
        }
        printf("\033[H");
        for (int i = 0; i < SCREEN_HEIGHT; ++i) {
            for (int j = 0; j < SCREEN_WIDTH; ++j) {
                putchar(output[i * SCREEN_WIDTH + j]);
            }
            putchar('\n');
        }
        A += 0.04;
        B += 0.02;
        sleep_ms(30);
#ifdef _WIN32
        if (_kbhit()) {
            char c = _getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#else
        if (kbhit()) {
            char c = getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#endif
    }
    printf("\nBuckyball animation stopped.\n");
}

// High-quality shaded triangular prism animation
void animate_triangular_prism() {
    float A = 0, B = 0;
    char output[SCREEN_WIDTH * SCREEN_HEIGHT];
    float zbuffer[SCREEN_WIDTH * SCREEN_HEIGHT];
    const char *lum = ".,-~:;=!*#$@";
    int stop = 0;
    printf("\033[2J");
    // Triangular prism: two triangles (top and bottom), three rectangles (sides)
    float top[3][3] = {{-1,1,0},{1,1,0},{0,1,sqrt(3)}};
    float bot[3][3] = {{-1,-1,0},{1,-1,0},{0,-1,sqrt(3)}};
    // Sides: connect corresponding vertices
    while (!stop) {
        memset(output, ' ', SCREEN_WIDTH * SCREEN_HEIGHT);
        memset(zbuffer, 0, sizeof(zbuffer));
        // Top and bottom triangles
        for (int t = 0; t < 2; ++t) {
            float (*tri)[3] = t ? bot : top;
            for (float u = 0; u <= 1.0; u += 0.03) {
                for (float v = 0; v <= 1.0 - u; v += 0.03) {
                    float w = 1.0 - u - v;
                    float x = w * tri[0][0] + u * tri[1][0] + v * tri[2][0];
                    float y = w * tri[0][1] + u * tri[1][1] + v * tri[2][1];
                    float z = w * tri[0][2] + u * tri[1][2] + v * tri[2][2];
                    // Normal for triangle
                    float ux = tri[1][0] - tri[0][0], uy = tri[1][1] - tri[0][1], uz = tri[1][2] - tri[0][2];
                    float vx_ = tri[2][0] - tri[0][0], vy_ = tri[2][1] - tri[0][1], vz_ = tri[2][2] - tri[0][2];
                    float nx = uy * vz_ - uz * vy_;
                    float ny = uz * vx_ - ux * vz_;
                    float nz = ux * vy_ - uy * vx_;
                    // Rotate point and normal
                    float x1 = x * cos(B) + z * sin(B);
                    float y1 = y;
                    float z1 = -x * sin(B) + z * cos(B);
                    float x2 = x1;
                    float y2 = y1 * cos(A) - z1 * sin(A);
                    float z2 = y1 * sin(A) + z1 * cos(A);
                    float nx1 = nx * cos(B) + nz * sin(B);
                    float ny1 = ny;
                    float nz1 = -nx * sin(B) + nz * cos(B);
                    float nx2 = nx1;
                    float ny2 = ny1 * cos(A) - nz1 * sin(A);
                    float nz2 = ny1 * sin(A) + nz1 * cos(A);
                    float K = 8;
                    float ooz = 1 / (K2 + z2);
                    int xp = (int)(SCREEN_WIDTH / 2 + K1 * ooz * x2);
                    int yp = (int)(SCREEN_HEIGHT / 2 - K1 * ooz * y2);
                    int idx = xp + yp * SCREEN_WIDTH;
                    float Lx = 0, Ly = 1, Lz = -1;
                    float norm = sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);
                    float L = (nx2 * Lx + ny2 * Ly + nz2 * Lz) / (norm * sqrt(Lx * Lx + Ly * Ly + Lz * Lz));
                    if (idx >= 0 && idx < SCREEN_WIDTH * SCREEN_HEIGHT && L > 0) {
                        if (ooz > zbuffer[idx]) {
                            zbuffer[idx] = ooz;
                            int lum_idx = (int)(L * 10);
                            if (lum_idx < 0) lum_idx = 0;
                            if (lum_idx > 11) lum_idx = 11;
                            output[idx] = lum[lum_idx];
                        }
                    }
                }
            }
        }
        // Three rectangular sides
        for (int s = 0; s < 3; ++s) {
            float *t0 = top[s], *t1 = top[(s+1)%3];
            float *b0 = bot[s], *b1 = bot[(s+1)%3];
            for (float u = 0; u <= 1.0; u += 0.03) {
                for (float v = 0; v <= 1.0; v += 0.03) {
                    float x = (1-u)*(1-v)*t0[0] + u*(1-v)*t1[0] + u*v*b1[0] + (1-u)*v*b0[0];
                    float y = (1-u)*(1-v)*t0[1] + u*(1-v)*t1[1] + u*v*b1[1] + (1-u)*v*b0[1];
                    float z = (1-u)*(1-v)*t0[2] + u*(1-v)*t1[2] + u*v*b1[2] + (1-u)*v*b0[2];
                    // Normal for side
                    float ux = t1[0] - t0[0], uy = t1[1] - t0[1], uz = t1[2] - t0[2];
                    float vx_ = b0[0] - t0[0], vy_ = b0[1] - t0[1], vz_ = b0[2] - t0[2];
                    float nx = uy * vz_ - uz * vy_;
                    float ny = uz * vx_ - ux * vz_;
                    float nz = ux * vy_ - uy * vx_;
                    // Rotate point and normal
                    float x1 = x * cos(B) + z * sin(B);
                    float y1 = y;
                    float z1 = -x * sin(B) + z * cos(B);
                    float x2 = x1;
                    float y2 = y1 * cos(A) - z1 * sin(A);
                    float z2 = y1 * sin(A) + z1 * cos(A);
                    float nx1 = nx * cos(B) + nz * sin(B);
                    float ny1 = ny;
                    float nz1 = -nx * sin(B) + nz * cos(B);
                    float nx2 = nx1;
                    float ny2 = ny1 * cos(A) - nz1 * sin(A);
                    float nz2 = ny1 * sin(A) + nz1 * cos(A);
                    float K = 8;
                    float ooz = 1 / (K2 + z2);
                    int xp = (int)(SCREEN_WIDTH / 2 + K1 * ooz * x2);
                    int yp = (int)(SCREEN_HEIGHT / 2 - K1 * ooz * y2);
                    int idx = xp + yp * SCREEN_WIDTH;
                    float Lx = 0, Ly = 1, Lz = -1;
                    float norm = sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);
                    float L = (nx2 * Lx + ny2 * Ly + nz2 * Lz) / (norm * sqrt(Lx * Lx + Ly * Ly + Lz * Lz));
                    if (idx >= 0 && idx < SCREEN_WIDTH * SCREEN_HEIGHT && L > 0) {
                        if (ooz > zbuffer[idx]) {
                            zbuffer[idx] = ooz;
                            int lum_idx = (int)(L * 10);
                            if (lum_idx < 0) lum_idx = 0;
                            if (lum_idx > 11) lum_idx = 11;
                            output[idx] = lum[lum_idx];
                        }
                    }
                }
            }
        }
        printf("\033[H");
        for (int i = 0; i < SCREEN_HEIGHT; ++i) {
            for (int j = 0; j < SCREEN_WIDTH; ++j) {
                putchar(output[i * SCREEN_WIDTH + j]);
            }
            putchar('\n');
        }
        A += 0.04;
        B += 0.02;
        sleep_ms(30);
#ifdef _WIN32
        if (_kbhit()) {
            char c = _getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#else
        if (kbhit()) {
            char c = getch();
            if (c == 'S' || c == 's') stop = 1;
        }
#endif
    }
    printf("\nTriangular prism animation stopped.\n");
}

int main() {
    int choice = 0;
    printf("\nASCII 3D Shape Renderer\n");
    printf("1. Donut\n2. Sphere\n3. Cube\n4. Pyramid\n5. Cone\n6. Cylinder\n7. Hourglass\n8. Icosahedron\n9. Buckyball\n10. Triangular Prism\nChoose a shape (1-10): ");
    scanf("%d", &choice);
    switch (choice) {
        case 1:
            animate_donut();
            break;
        case 2:
            animate_sphere();
            break;
        case 3:
            animate_cube();
            break;
        case 4:
            animate_pyramid();
            break;
        case 5:
            animate_cone();
            break;
        case 6:
            animate_cylinder();
            break;
        case 7:
            animate_hourglass();
            break;
        case 8:
            animate_icosahedron();
            break;
        case 9:
            animate_buckyball();
            break;
        case 10:
            animate_triangular_prism();
            break;
        default:
            printf("Invalid choice. Exiting.\n");
    }
    return 0;
} 