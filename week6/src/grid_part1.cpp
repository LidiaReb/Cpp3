#include<iostream>

template <typename T>
class Grid final {
public:
    using value_type = T;
    using size_type = unsigned;

private:
    T * const data;
    size_type const y_size, x_size;

public:
    Grid(T *data , size_type y_size , size_type x_size): data(data), y_size(y_size), x_size(x_size) {}

    //конструктор (неявное преобразование типов T -> Grid<T>)
    Grid(T const &t): Grid(y_size, x_size, t) {} 
    
    //конструктор со значением default T по всей сетке
    Grid(size_type y_size, size_type x_size): Grid(y_size, x_size, T()) {}
    
    //конструктор со значением t по всей сетке
    Grid(size_type y_size, size_type x_size, T const &t):
        data(reinterpret_cast<T*>( operator new(sizeof(T) * y_size * x_size) )), 
        y_size(y_size), x_size(x_size) {
            for (size_type idx = 0u; idx <  y_size*x_size; ++idx) {
                new ((char*)data + idx * sizeof(T)) T(t);
            }     
    }

    //конструктор копирования 
    Grid(Grid<T> const &copy):
        data(new T[copy.y_size * copy.x_size]), 
        y_size(copy.y_size), 
        x_size(copy.x_size) {
        for (std::size_t idx = 0u; idx != y_size*x_size; ++idx) {
            data[idx] = copy.data[idx];
        }    
    }

    //конструктор перемещения
    Grid(Grid<T> &&copy):
        data(copy.data), 
        y_size(copy.y_size), 
        x_size(copy.x_size) {
        copy.data = nullptr;
        copy.y_size = 0;
        copy.x_size = 0; 
    }

    //оператор присваивания копированием
    Grid<T>& operator=(Grid<T> const &src) {
        for (std::size_t idx = 0u; idx != y_size*x_size; ++idx) {
            data[idx] = src[idx];
        } 
        return *this;
    }

    //оператор присваивания перемещением
    Grid<T>& operator=(Grid<T> &&src) {  
        for (std::size_t idx = 0u; idx != y_size*x_size; ++idx) {
            data[idx] = src[idx];
        } 
        src.data = nullptr;
        src.y_size = 0;
        src.x_size = 0; 
        return *this;
    }

    T operator()(size_type y_idx, size_type x_idx) const {
        return data[y_idx * x_size + x_idx];
    }

    T& operator()(size_type y_idx, size_type x_idx){
        return data[y_idx * x_size + x_idx];
    }

    //оператор присваивания копированием одинакового элемента во все ячейки сетки
    Grid<T>& operator=(T const &t){
        for (auto it = data, end = data + x_size * y_size; it != end; ++it) 
            *it = t;
        return *this;
    }

    size_type get_y_size() const { return y_size; }
    size_type get_x_size() const { return x_size; }

    void print() const {
        for (size_type i = 0; i < y_size; ++i) {
            std::cout << i << ":\t"; 
            for (size_type j = 0; j < x_size; ++j) {
                std::cout << data[i * x_size + j] << '\t';
            }
            std::cout << '\n';
        }
    }

    ~Grid() {
    for(size_type idx = 0; idx < y_size * x_size; ++idx) {
            data[idx].~T();
        }
        delete[] data;
    } 
};



