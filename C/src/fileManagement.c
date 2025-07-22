#include <sys/stat.h>  // For mkdir() and stat()
#include <sys/types.h> // For mode_t
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#define BUFFER_SIZE 4096  // Size of the buffer used to copy

#include "../include/fileManagement.h"
#include "../include/sysVar.h"

int createFolder(const char* dir) {

    struct stat st = {0};

    // Check if the directory already exists
    if (stat(dir, &st) == -1) {

        // Create the directory
        if (mkdir(dir, 0700) == 0) {
            printf("Directory '%s' created successfully.\n", dir);
        
        } else {
            perror("Error creating directory");
            return 1;
        }
    } else {
        printf("Directory '%s' already exists.\n", dir);

    }

    return 0;
}


int createDirectories(const char *path) {
    
    int check = 0;
    char temp[256];

    strcpy(temp, path);

    // Pointing to the address of the first element of temp
    char* pos = temp;

    while ((pos = strchr(pos, '/')) != NULL) {
        
        *pos = '\0';  // Temporarily terminate the string at this position

        // Check if this part of the path exists
        check = createFolder(temp);
        if (check == 1){ return 1; }

        *pos = '/';  // Restore the character and move on
        pos++;
    }

    // Create the final directory or subdirectory
    check = createFolder(temp);
    if (check == 1){ return 1; }

    return 0;
}


void float2stringName(double number, char* buffer, size_t size) {

    // Step 1: Format the float to a string with 6 decimal places
    snprintf(buffer, size, "%.8f", number);

    // Step 2: Replace the decimal point with 'p'
    // char *decimal_point = strchr(buffer, '.');

    // if (decimal_point != NULL) {
    //     *decimal_point = 'p';
    // }

    return;
}


void createBaseName(const double Tstar, const double phi, char* baseName, size_t size) {

    char Tstar_str[20];
    char phi_str[20];

    float2stringName(Tstar, Tstar_str, sizeof(Tstar_str));
    float2stringName(phi, phi_str, sizeof(phi_str));

    // Use snprintf to format and concatenate the strings
    snprintf(baseName, size, "T_%s_phi_%s", Tstar_str, phi_str);

    return;
}

void createBaseFolder(const double Tstar, const double phi, char* baseFolder, size_t size) {

    char Tstar_str[20];
    char phi_str[20];

    float2stringName(Tstar, Tstar_str, sizeof(Tstar_str));
    float2stringName(phi, phi_str, sizeof(phi_str));

    // Use snprintf to format and concatenate the strings
    snprintf(baseFolder, size, "T_%s/phi_%s", Tstar_str, phi_str);

    return;
}

void createPath(const char* folder, const char* baseName, const char* prefix, const char* extension, char* path, size_t size) {

    char filename[200];

    snprintf(filename, sizeof(filename), "%s_%s.%s", prefix, baseName, extension);

    snprintf(path, size, "%s/%s", folder, filename);

    return;
}

void writeTensor(char* filename, double* x, tensor3d* A1, SysVar sysVar) {

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error: Could not open hr matrix output file :'( . \n");
        return;
        // return EXIT_FAILURE;
    }
    
    size_t N1 = A1->mats_number;
    size_t N2 = A1->mats_row_size;
    size_t N3 = A1->mats_col_size;

    // Write the matrix to the file
    for (size_t o=0; o<N1; o++) {
        fprintf(file, "%.10e ", x[o]);

        for (size_t i=0; i<N2; i++) {
            for (size_t j=0; j<N3; j++) {

                fprintf(file, "%.10e ", A1->tensor[o][i][j]);
            }
        }

        fprintf(file, "\n");
    }

    fclose(file);

    return;
}

void writeTensorFiles(const char* baseFolder, const char* baseName, const char* inputVarFile, 
                      MainVar mainVar, SysVar sysVar) {

    char path[300];

    createDirectories(baseFolder);

    createPath(baseFolder, baseName, "hrTensor", "txt", path, sizeof(path));
    writeTensor(path, mainVar.x, mainVar.hrMat, sysVar);

    createPath(baseFolder, baseName, "crTensor", "txt", path, sizeof(path));
    writeTensor(path, mainVar.x, mainVar.crMat, sysVar);

    createPath(baseFolder, baseName, "SkTensor", "txt", path, sizeof(path));
    writeTensor(path, mainVar.x, mainVar.SkMat, sysVar);

    createPath(baseFolder, baseName, "CkTensor", "txt", path, sizeof(path));
    writeTensor(path, mainVar.x, mainVar.CkMat, sysVar);

    createPath(baseFolder, baseName, inputVarFile, "txt", path, sizeof(path));
    int checkCopy = copyFile(inputVarFile, path);

    if (checkCopy == 0) {
        printf("File copied successfully!\n");
    } else {
        printf("Failed to copy the file.\n");
    }

    return;
}

void printdPdRho(double dPdRho) {

    printf("\n----------------------------\n");
    printf("dPdRho = %.10lf \n", dPdRho);
    printf("----------------------------\n");

    return;
}


int copyFile(const char* sourcePath, const char* destinationPath) {
    FILE *sourceFile, *destFile;
    char buffer[BUFFER_SIZE];
    size_t bytesRead;

    // Open the source file in read mode (binary mode for binary files)
    sourceFile = fopen(sourcePath, "rb");
    if (sourceFile == NULL) {
        perror("Error opening source file");
        return 1;
    }

    // Open the destination file in write mode (binary mode for binary files)
    destFile = fopen(destinationPath, "wb");
    if (destFile == NULL) {
        perror("Error opening destination file");
        fclose(sourceFile);  // Close the source file before returning
        return 1;
    }

    // Copy data from source to destination in chunks
    while ((bytesRead = fread(buffer, 1, BUFFER_SIZE, sourceFile)) > 0) {
        fwrite(buffer, 1, bytesRead, destFile);
    }

    // Close both files
    fclose(sourceFile);
    fclose(destFile);

    return 0;
}

char* concatenate(const char* folder, const char* name) {
    
    // Calculate the needed size (+1 for the possible '/' and +1 for the null terminator)
    size_t folder_len = strlen(folder);
    size_t name_len = strlen(name);
    size_t total_len = folder_len + name_len + 2;

    // Allocate memory for the result
    char* fullPath = (char*)malloc(total_len);
    if (fullPath == NULL) {
        perror("Unable to allocate memory");
        return NULL;
    }

    // Copy the folder path
    strcpy(fullPath, folder);

    // Add a separator '/' if needed
    if (folder[folder_len - 1] != '/' && folder[folder_len - 1] != '\\') {
        strcat(fullPath, "/");
    }

    // Append the name
    strcat(fullPath, name);

    return fullPath;
}

const char* getFileName(const char* path) {
    const char* slash = strrchr(path, '/');
    const char* backslash = strrchr(path, '\\');
    if (slash || backslash) {
        return slash > backslash ? slash + 1 : backslash + 1;
    }
    return path; // If no slash or backslash, the path is already the file name
}