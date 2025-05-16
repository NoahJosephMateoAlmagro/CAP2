//==================================================================================================
// Written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is distributed
// without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication along
// with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==================================================================================================

#include <float.h>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <sstream>
#include <fstream>
#include <omp.h>

#include "Camera.h"
#include "Object.h"
#include "Scene.h"
#include "Sphere.h"
#include "Diffuse.h"
#include "Metallic.h"
#include "Crystalline.h"

#include "random.h"
#include "utils.h"
#include <chrono> 


Scene loadObjectsFromFile(const std::string& filename) {
	std::ifstream file(filename);
	std::string line;

	Scene list;

	if (file.is_open()) {
		while (std::getline(file, line)) {
			std::stringstream ss(line);
			std::string token;
			std::vector<std::string> tokens;

			while (ss >> token) {
				tokens.push_back(token);
			}

			if (tokens.empty()) continue; // Línea vacía

			// Esperamos al menos la palabra clave "Object"
			if (tokens[0] == "Object" && tokens.size() >= 12) { // Mínimo para Sphere y un material con 1 float
				// Parsear la esfera
				if (tokens[1] == "Sphere" && tokens[2] == "(" && tokens[7] == ")") {
					try {
						float sx = std::stof(tokens[3].substr(tokens[3].find('(') + 1, tokens[3].find(',') - tokens[3].find('(') - 1));
						float sy = std::stof(tokens[4].substr(0, tokens[4].find(',')));
						float sz = std::stof(tokens[5].substr(0, tokens[5].find(',')));
						float sr = std::stof(tokens[6]);

						// Parsear el material del último objeto creado

						if (tokens[8] == "Crystalline" && tokens[9] == "(" && tokens[11].back() == ')') {
							float ma = std::stof(tokens[10]);
							list.add(new Object(
								new Sphere(Vec3(sx, sy, sz), sr),
								new Crystalline(ma)
							));
							std::cout << "Crystaline" << sx << " " << sy << " " << sz << " " << sr << " " << ma << "\n";
						}
						else if (tokens[8] == "Metallic" && tokens.size() == 15 && tokens[9] == "(" && tokens[14] == ")") {
							float ma = std::stof(tokens[10].substr(tokens[10].find('(') + 1, tokens[10].find(',') - tokens[10].find('(') - 1));
							float mb = std::stof(tokens[11].substr(0, tokens[11].find(',')));
							float mc = std::stof(tokens[12].substr(0, tokens[12].find(',')));
							float mf = std::stof(tokens[13].substr(0, tokens[13].length() - 1));
							list.add(new Object(
								new Sphere(Vec3(sx, sy, sz), sr),
								new Metallic(Vec3(ma, mb, mc), mf)
							));
							std::cout << "Metallic" << sx << " " << sy << " " << sz << " " << sr << " " << ma << " " << mb << " " << mc << " " << mf << "\n";
						}
						else if (tokens[8] == "Diffuse" && tokens.size() == 14 && tokens[9] == "(" && tokens[13].back() == ')') {
							float ma = std::stof(tokens[10].substr(tokens[10].find('(') + 1, tokens[10].find(',') - tokens[10].find('(') - 1));
							float mb = std::stof(tokens[11].substr(0, tokens[11].find(',')));
							float mc = std::stof(tokens[12].substr(0, tokens[12].find(',')));
							list.add(new Object(
								new Sphere(Vec3(sx, sy, sz), sr),
								new Diffuse(Vec3(ma, mb, mc))
							));
							std::cout << "Diffuse" << sx << " " << sy << " " << sz << " " << sr << " " << ma << " " << mb << " " << mc << "\n";
						}
						else {
							std::cerr << "Error: Material desconocido o formato incorrecto en la línea: " << line << std::endl;
						}
					}
					catch (const std::invalid_argument& e) {
						std::cerr << "Error: Conversión inválida en la línea: " << line << " - " << e.what() << std::endl;
					}
					catch (const std::out_of_range& e) {
						std::cerr << "Error: Valor fuera de rango en la línea: " << line << " - " << e.what() << std::endl;
					}
				}
				else {
					std::cerr << "Error: Formato de esfera incorrecto en la línea: " << line << std::endl;
				}
			}
			else {
				std::cerr << "Error: Formato de objeto incorrecto en la línea: " << line << std::endl;
			}
		}
		file.close();
	}
	else {
		std::cerr << "Error: No se pudo abrir el archivo: " << filename << std::endl;
	}
	return list;
}

Scene randomScene()
{
	int n = 500;
	const int RANGO = 11; //Añadido para modificar el número de esferas.
	Scene list;
	list.add(new Object(
		new Sphere(Vec3(0, -1000, 0), 1000),
		new Diffuse(Vec3(0.5, 0.5, 0.5))));

	for (int a = -RANGO; a < RANGO; a++)
	{
		for (int b = -RANGO; b < RANGO; b++)
		{
			float choose_mat = Mirandom();
			Vec3 center(a + 0.9f * Mirandom(), 0.2f, b + 0.9f * Mirandom());
			if ((center - Vec3(4, 0.2f, 0)).length() > 0.9f)
			{
				if (choose_mat < 0.8f)
				{ // diffuse
					list.add(new Object(
						new Sphere(center, 0.2f),
						new Diffuse(Vec3(Mirandom() * Mirandom(),
										 Mirandom() * Mirandom(),
										 Mirandom() * Mirandom()))));
				}
				else if (choose_mat < 0.95f)
				{ // metal
					list.add(new Object(
						new Sphere(center, 0.2f),
						new Metallic(Vec3(0.5f * (1 + Mirandom()),
										  0.5f * (1 + Mirandom()),
										  0.5f * (1 + Mirandom())),
									 0.5f * Mirandom())));
				}
				else
				{ // glass
					list.add(new Object(
						new Sphere(center, 0.2f),
						new Crystalline(1.5f)));
				}
			}
		}
	}

	list.add(new Object(
		new Sphere(Vec3(0, 1, 0), 1.0),
		new Crystalline(1.5f)));
	list.add(new Object(
		new Sphere(Vec3(-4, 1, 0), 1.0f),
		new Diffuse(Vec3(0.4f, 0.2f, 0.1f))));
	list.add(new Object(
		new Sphere(Vec3(4, 1, 0), 1.0f),
		new Metallic(Vec3(0.7f, 0.6f, 0.5f), 0.0f)));

	return list;
}

void rayTracingCPU_OMP_Bloque(unsigned char *img, int w, int h, int ns = 10, int px = 0, int py = 0, int pw = -1, int ph = -1)
{
    if (pw == -1) pw = w;
    if (ph == -1) ph = h;

    std::cout << "[DEBUG] Iniciando generación de escena...\n";
    Scene world = randomScene();
    std::cout << "[DEBUG] Escena generada.\n";

    world.setSkyColor(Vec3(0.5f, 0.7f, 1.0f));
    world.setInfColor(Vec3(1.0f, 1.0f, 1.0f));

    Vec3 lookfrom(13, 2, 3);
    Vec3 lookat(0, 0, 0);
    float dist_to_focus = 10.0;
    float aperture = 0.1f;

    Camera cam(lookfrom, lookat, Vec3(0, 1, 0), 20, float(w) / float(h), aperture, dist_to_focus); 

    std::cout << "[DEBUG] Comenzando render...\n";

	const int BLOCK_SIZE = 16;

    int numBlocksX = (w + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int numBlocksY = (h + BLOCK_SIZE - 1) / BLOCK_SIZE;

    #pragma omp parallel for collapse(2) schedule(dynamic)
	
    for (int blockY = 0; blockY < numBlocksY; blockY++) {
        for (int blockX = 0; blockX < numBlocksX; blockX++) {

            int thread_id = omp_get_thread_num();
            printf("[THREAD %d] Renderizando bloque (%d, %d)\n", thread_id, blockX, blockY);

            int x_start = blockX * BLOCK_SIZE;
            int y_start = blockY * BLOCK_SIZE;

            int x_end = std::min(x_start + BLOCK_SIZE, w);
            int y_end = std::min(y_start + BLOCK_SIZE, h);

            for (int j = y_start; j < y_end; j++) {
                for (int i = x_start; i < x_end; i++) {

                    Vec3 col(0, 0, 0);
                    for (int s = 0; s < ns; s++) {
                        float u = float(i + Mirandom()) / float(w);
                        float v = float(j + Mirandom()) / float(h);
                        Ray r = cam.get_ray(u, v);
                        col += world.getSceneColor(r);
                    }

                    col /= float(ns);
                    col = Vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));

                    int idx = j * w * 3 + i * 3;
                    img[idx + 2] = char(255.99 * col[0]);
                    img[idx + 1] = char(255.99 * col[1]);
                    img[idx + 0] = char(255.99 * col[2]);
                }
            }
        }
    }

}

int main(int argc, char *argv[])
{
    int w = 1200;
    int h = 1200;
    int ns = 10;

    int size = sizeof(unsigned char) * w * h * 3;
    unsigned char *data = (unsigned char *)calloc(size, 1);

    int max_threads = omp_get_max_threads();
    printf("[DEBUG] Usando OpenMP con hasta %d hilos.\n", max_threads);

    double start_time = omp_get_wtime();

    rayTracingCPU_OMP_Bloque(data, w, h, ns);

    double end_time = omp_get_wtime();

    writeBMP("imgOMPBloque.bmp", data, w, h);
    printf("[DEBUG] Imagen creada: imgOMPFilas.bmp\n");
    printf("Tiempo de ejecución con OpenMP (por filas): %f segundos.\n", end_time - start_time);

    free(data);
    return 0;
}