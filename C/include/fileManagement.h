#ifndef FILEMANAGEMENT_H
#define FILEMANAGEMENT_H

#include "sysVar.h"

int createFolder(const char* dir);
int createDirectories(const char *path);
void float2stringName(double number, char* buffer, size_t size);
void createBaseName(const double Tstar, const double phi, char* baseName, size_t size);
void createBaseFolder(const double Tstar, const double phi, char* baseFolder, size_t size);
void createPath(const char* folder, const char* baseName, const char* prefix, const char* extension, char* path, size_t size);
void writeTensor(char* filename, double* x, tensor3d* A1, SysVar sysVar);
void writeTensorFiles(const char* baseFolder, const char* baseName, const char* inputVarFile, 
                      MainVar mainVar, SysVar sysVar);
void printdPdRho(double dPdRho);
int copyFile(const char* sourcePath, const char* destinationPath);

char* concatenate(const char* folder, const char* name);
const char* getFileName(const char* path);

#endif // FILEMANAGEMENT_H