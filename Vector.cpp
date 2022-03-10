#include <iostream>
#include "Vector.h"

Vector::Vector(double x, double y, double z, std::string name, int type)
{            
    this->x = x;
    this->y = y;
    this->z = z;
    this->name = name;
    this->type = type;
}
void Vector::get_coord() 
{
    std::cout << this->name << ' ' << this->x << ' ' << this->y << ' ' << this->z << std::endl;
}

double Vector::getX() 
{
    return this->x;
}

double Vector::getY()
{
    return this->y;
}

double Vector::getZ() 
{
    return this->z;
}

int Vector::getType()
{
    return this->type;
}

std::string Vector::getName()
{
    return this->name;
}