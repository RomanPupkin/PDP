#include <iostream>

class Vector
{
    private:
        int type;
        double x, y, z;
        std::string name;
    public:
        double getX();
        double getY();
        double getZ();
        int getType();
        std::string getName();
        Vector& operator*(double coeff);
        Vector(double x, double y, double z, std::string name, int type);
        ~Vector() = default;
        void get_coord();
};
