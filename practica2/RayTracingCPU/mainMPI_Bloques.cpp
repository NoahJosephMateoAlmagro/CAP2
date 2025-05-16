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
#include <iomanip>

#include "Camera.h"
#include "Object.h"
#include "Scene.h"
#include "Sphere.h"
#include "Diffuse.h"
#include "Metallic.h"
#include "Crystalline.h"

#include "random.h"
#include "utils.h"

#include <mpi.h>
#include <omp.h>

Scene loadObjectsFromFile(const std::string &filename)
{
	std::ifstream file(filename);
	std::string line;

	Scene list;

	if (file.is_open())
	{
		while (std::getline(file, line))
		{
			std::stringstream ss(line);
			std::string token;
			std::vector<std::string> tokens;

			while (ss >> token)
			{
				tokens.push_back(token);
			}

			if (tokens.empty())
				continue; // L�nea vac�a

			// Esperamos al menos la palabra clave "Object"
			if (tokens[0] == "Object" && tokens.size() >= 12)
			{ // M�nimo para Sphere y un material con 1 float
				// Parsear la esfera
				if (tokens[1] == "Sphere" && tokens[2] == "(" && tokens[7] == ")")
				{
					try
					{
						float sx = std::stof(tokens[3].substr(tokens[3].find('(') + 1, tokens[3].find(',') - tokens[3].find('(') - 1));
						float sy = std::stof(tokens[4].substr(0, tokens[4].find(',')));
						float sz = std::stof(tokens[5].substr(0, tokens[5].find(',')));
						float sr = std::stof(tokens[6]);

						// Parsear el material del �ltimo objeto creado

						if (tokens[8] == "Crystalline" && tokens[9] == "(" && tokens[11].back() == ')')
						{
							float ma = std::stof(tokens[10]);
							list.add(new Object(
								new Sphere(Vec3(sx, sy, sz), sr),
								new Crystalline(ma)));
							//std::cout << "Crystaline" << sx << " " << sy << " " << sz << " " << sr << " " << ma << "\n";
						}
						else if (tokens[8] == "Metallic" && tokens.size() == 15 && tokens[9] == "(" && tokens[14] == ")")
						{
							float ma = std::stof(tokens[10].substr(tokens[10].find('(') + 1, tokens[10].find(',') - tokens[10].find('(') - 1));
							float mb = std::stof(tokens[11].substr(0, tokens[11].find(',')));
							float mc = std::stof(tokens[12].substr(0, tokens[12].find(',')));
							float mf = std::stof(tokens[13].substr(0, tokens[13].length() - 1));
							list.add(new Object(
								new Sphere(Vec3(sx, sy, sz), sr),
								new Metallic(Vec3(ma, mb, mc), mf)));
							//std::cout << "Metallic" << sx << " " << sy << " " << sz << " " << sr << " " << ma << " " << mb << " " << mc << " " << mf << "\n";
						}
						else if (tokens[8] == "Diffuse" && tokens.size() == 14 && tokens[9] == "(" && tokens[13].back() == ')')
						{
							float ma = std::stof(tokens[10].substr(tokens[10].find('(') + 1, tokens[10].find(',') - tokens[10].find('(') - 1));
							float mb = std::stof(tokens[11].substr(0, tokens[11].find(',')));
							float mc = std::stof(tokens[12].substr(0, tokens[12].find(',')));
							list.add(new Object(
								new Sphere(Vec3(sx, sy, sz), sr),
								new Diffuse(Vec3(ma, mb, mc))));
							//std::cout << "Diffuse" << sx << " " << sy << " " << sz << " " << sr << " " << ma << " " << mb << " " << mc << "\n";
						}
						else
						{
							std::cerr << "Error: Material desconocido o formato incorrecto en la l�nea: " << line << std::endl;
						}
					}
					catch (const std::invalid_argument &e)
					{
						std::cerr << "Error: Conversi�n inv�lida en la l�nea: " << line << " - " << e.what() << std::endl;
					}
					catch (const std::out_of_range &e)
					{
						std::cerr << "Error: Valor fuera de rango en la l�nea: " << line << " - " << e.what() << std::endl;
					}
				}
				else
				{
					std::cerr << "Error: Formato de esfera incorrecto en la l�nea: " << line << std::endl;
				}
			}
			else
			{
				std::cerr << "Error: Formato de objeto incorrecto en la l�nea: " << line << std::endl;
			}
		}
		file.close();
	}
	else
	{
		std::cerr << "Error: No se pudo abrir el archivo: " << filename << std::endl;
	}
	return list;
}

Scene randomScene()
{
	int n = 500;
	const int RANGO = 11; // Añadido para modificar el número de esferas.
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

// Función para guardar la escena en un txt (fuerza a añadir getters y setters en los materiales)
void saveSceneToFile(const std::string &filename, const Scene &scene)
{
	std::ofstream file(filename);
	if (!file.is_open())
	{
		std::cerr << "Error: No se pudo abrir el archivo para guardar: " << filename << std::endl;
		return;
	}

	file << std::fixed << std::setprecision(1); // Para guardar los enteros como X.0

	for (const auto &obj : scene.getObjects())
	{
		Sphere *sphere = dynamic_cast<Sphere *>(obj->getShape());
		if (!sphere)
			continue;

		Vec3 center = sphere->getCenter();
		float radius = sphere->getRadius();

		Material *mat = obj->getMaterial();
		if (auto *diffuse = dynamic_cast<Diffuse *>(mat))
		{
			Vec3 color = diffuse->getColor();
			file << "Object Sphere ( ("
				 << center[0] << ", " << center[1] << ", " << center[2]
				 << "), " << radius << " ) Diffuse ( ("
				 << color[0] << ", " << color[1] << ", " << color[2] << ") )\n";
		}
		else if (auto *m = dynamic_cast<Metallic *>(mat))
		{
			Vec3 albedo = m->getAlbedo();
			float fuzz = m->getFuzziness();
			file << "Object Sphere ( ("
				 << center[0] << ", " << center[1] << ", " << center[2]
				 << "), " << radius << " ) Metallic ( ("
				 << albedo[0] << ", " << albedo[1] << ", " << albedo[2]
				 << "), " << fuzz << " )\n";
		}
		else if (auto *c = dynamic_cast<Crystalline *>(mat))
		{
			float ref_idx = c->getRefIdx();
			file << "Object Sphere ( ("
				 << center[0] << ", " << center[1] << ", " << center[2]
				 << "), " << radius << " ) Crystalline ( "
				 << ref_idx << " )\n";
		}
		else
		{
			std::cerr << "Error: Material desconocido al guardar.\n";
		}
	}

	file.close();
}

void rayTracingCPU_MPI(unsigned char *img, int w, int h, int ns = 10, int px = 0, int py = 0, int pw = -1, int ph = -1)
{
	if (pw == -1)
		pw = w;
	if (ph == -1)
		ph = h;
	int patch_w = pw - px;

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	Scene world;

	if (rank == 0)
	{

		world = randomScene();
		// Scene world = loadObjectsFromFile("Scene1.txt");

		saveSceneToFile("Scene2.txt", world);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	world = loadObjectsFromFile("Scene2.txt"); //El proceso 0 lo vuelve a leer para evitar inconsistencia
	world.setSkyColor(Vec3(0.5f, 0.7f, 1.0f));
	world.setInfColor(Vec3(1.0f, 1.0f, 1.0f));

	Vec3 lookfrom(13, 2, 3);
	Vec3 lookat(0, 0, 0);
	float dist_to_focus = 10.0;
	float aperture = 0.1f;

	Camera cam(lookfrom, lookat, Vec3(0, 1, 0), 20, float(w) / float(h), aperture, dist_to_focus);

	// Elegir la mejor factorización (bx * by = size)
	int bx = 1, by = size;
	for (int i = 1; i <= sqrt(size); ++i)
	{
		if (size % i == 0)
		{
			bx = i;
			by = size / i;
		}
	}

	int proc_x = rank % bx;
	int proc_y = rank / bx;

	int base_w = w / bx, extra_w = w % bx;
	int base_h = h / by, extra_h = h % by;

	int x_start = proc_x * base_w + std::min(proc_x, extra_w);
	int y_start = proc_y * base_h + std::min(proc_y, extra_h);

	int block_w = base_w + (proc_x < extra_w ? 1 : 0);
	int block_h = base_h + (proc_y < extra_h ? 1 : 0);

	printf("Proceso %d: x = %d, y = %d, ancho = %d, alto = %d\n",
		   rank, x_start, y_start, block_w, block_h);




	// Reservar imagen local
	unsigned char *local_img = (unsigned char *)calloc(block_w * block_h * 3, sizeof(unsigned char));

	// Renderizar la parte correspondiente
	for (int j = 0; j < block_h; ++j) {
        for (int i = 0; i < block_w; ++i) {
            int global_x = x_start + i;
            int global_y = y_start + j;

            Vec3 col(0, 0, 0);
            for (int s = 0; s < ns; ++s) {
                float u = float(global_x + Mirandom()) / float(w);
                float v = float(global_y + Mirandom()) / float(h);
                Ray r = cam.get_ray(u, v);
                col += world.getSceneColor(r);
            }

            col /= float(ns);
            col = Vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));

            int idx = (j * block_w + i) * 3;
            local_img[idx + 0] = char(255.99 * col[2]);
            local_img[idx + 1] = char(255.99 * col[1]);
            local_img[idx + 2] = char(255.99 * col[0]);
        }
    }


	// Preparar buffers en el proceso raíz

	if (rank == 0) {
        // Copia su propio bloque
        for (int j = 0; j < block_h; ++j) {
            memcpy(&img[((y_start + j) * w + x_start) * 3],
                   &local_img[j * block_w * 3],
                   block_w * 3);
        }

        // Recibir y ensamblar bloques de otros procesos
        for (int p = 1; p < size; ++p) {
            int px = p % bx;
            int py = p / bx;

            int xs = px * base_w + std::min(px, extra_w);
            int ys = py * base_h + std::min(py, extra_h);
            int bw = base_w + (px < extra_w ? 1 : 0);
            int bh = base_h + (py < extra_h ? 1 : 0);

            unsigned char * recv_buf = (unsigned char * )malloc(bw * bh * 3);
            MPI_Recv(recv_buf, bw * bh * 3, MPI_UNSIGNED_CHAR, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (int j = 0; j < bh; ++j) {
                memcpy(&img[((ys + j) * w + xs) * 3],
                       &recv_buf[j * bw * 3],
                       bw * 3);
            }

            free(recv_buf);
        }

    } else {
        MPI_Send(local_img, block_w * block_h * 3, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
    } 

	printf("Proceso %d: ha enviado %d bytes (pixeles: %d x %d)\n",
		   rank, block_w * block_h * 3, block_w, block_h);

	// Liberar memoria
	free(local_img);
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	int w = 160; // Ancho de la imagen
	int h = 160; // Alto de la imagen
	int ns = 10;  // Número de muestras por pixel

	int size = sizeof(unsigned char) * w * h * 3;
	unsigned char *data = (unsigned char *)calloc(size, 1);

	// Medir el tiempo de ejecución
	double start_time = MPI_Wtime(); // Iniciar el temporizador

	// Llamada a la función de ray tracing con MPI
	rayTracingCPU_MPI(data, w, h, ns);

	// Solo el proceso raíz escribe el archivo BMP
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		writeBMP("imgMPIImgBlocks.bmp", data, w, h);
		printf("Imagen creada.\n");
	}

	// Medir el tiempo después de la ejecución
	double end_time = MPI_Wtime(); // Finalizar el temporizador

	// Mostrar el tiempo de ejecución
	if (rank == 0)
	{
		printf("Tiempo de ejecución con SOLO MPI: %f segundos.\n", end_time - start_time);
	}

	free(data);
	MPI_Finalize();
	return 0;
}