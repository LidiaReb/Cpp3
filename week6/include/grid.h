#include <iostream>

template<typename T, std::size_t dim>
class Grid {
public:
    using value_type = T;
    using size_type = unsigned;
private:
    Grid<T, dim - 1>* buffer;
    size_type buffer_size;
public:
    //конструктор со значением default T или t по всей сетке
    template <typename ...Dims>
    Grid(size_type buffer_size, Dims... other); 

    //конструктор (неявное преобразование типов T -> Grid<T, dim> 1x...x1)
    Grid(T const &t);

    //конструктор по умолчанию (создает пустую сетку)
    Grid();

    //конструктор копирования
    Grid(Grid<T, dim> const &copy);

    //конструктор перемещения
    Grid(Grid<T, dim> &&copy); 

    //оператор присваивания копированием
    Grid<T, dim>& operator=(Grid<T, dim> const &src);

    //оператор присваивания перемещением
    Grid<T, dim>& operator=(Grid<T, dim> &&src);

    //деструктор
    ~Grid(); 

    //метод возвращает размер сетки по данному направлению
    size_type get_dim_size(std::size_t dim_num) const;

    //оператор индексирования
    Grid<T, dim - 1>& operator[] (size_type idx) const; 

    //другой оператор индексирования 
    template <typename ...Dims>
    T operator()(size_type idx, Dims... other) const;

    //другой оператор индексирования 
    template <typename ...Dims>
    T& operator()(size_type idx, Dims... other);

    //печать сетки
    void print() const ; 
};

template <typename T>
class Grid<T, 1> {
public:
    using value_type = T;
    using size_type = unsigned;
private:
    T* buffer;
    size_type buffer_size;
public:
    //конструктор со значением t по всей сетке
    Grid(std::size_t buffer_size, T const &t);

    //конструктор со значением default T по всей сетке
    Grid(std::size_t buffer_size);

    //конструктор (неявное преобразование типов T -> Grid<T, 1>)
    Grid(T const &t);

    //конструктор по умолчанию (создает пустую сетку)
    Grid();

    //конструктор копирования
    Grid(Grid<T, 1> const &copy);

    //конструктор перемещения
    Grid(Grid<T, 1> &&copy);

    //оператор присваивания копированием
    Grid<T, 1>& operator=(Grid<T, 1> const &src);

    //оператор присваивания перемещением
    Grid<T, 1>& operator=(Grid<T, 1> &&src);

    //деструктор
    ~Grid();

    //метод возвращает размер сетки по последнему направлению
    size_type get_dim_size(std::size_t dim_num) const;

    //оператор индексирования константный
    T& operator[] (size_type idx) const ;

    //оператор индексирования не константный
    T& operator[] (size_type idx);

    //другой оператор индексирования 
    T operator()(size_type idx) const ;

    //другой оператор индексирования 
    T& operator()(size_type idx);

    //печать сетки
    void print() const ;
};

template<typename T, std::size_t dim>
template<typename ...Dims>
Grid<T, dim>::Grid(size_type buffer_size, Dims... other) : 
    //конструктор со значением default T или t по всей сетке

    buffer(reinterpret_cast<Grid<T, dim - 1>*>( operator new(sizeof(Grid<T, dim - 1>) * buffer_size) )), 
    buffer_size(buffer_size) {
    for (size_type idx = 0; idx != buffer_size; ++idx)
        new ((char*)buffer + idx * sizeof(Grid<T, dim - 1>)) Grid<T, dim - 1>(other...);
}

template <typename T>
Grid<T, 1>::Grid(std::size_t buffer_size, T const &t) : 
    //конструктор со значением t по всей сетке

    buffer(reinterpret_cast<T*>( operator new(sizeof(T) * buffer_size) )), 
    buffer_size(buffer_size) {
    for (std::size_t idx = 0u; idx != buffer_size; ++idx) {
        new ((char*)buffer + idx * sizeof(T)) T(t);
    }
}

template <typename T>
Grid<T, 1>::Grid(std::size_t buffer_size) : 
    //конструктор со значением default T по всей сетке

    Grid<T, 1>(buffer_size, T()) { }

template <typename T, std::size_t dim>
Grid<T, dim>::Grid(T const &t) :
    //конструктор (неявное преобразование типов T -> Grid<T, dim> 1x...x1)

    buffer(reinterpret_cast<Grid<T, dim - 1>*>( operator new(sizeof(Grid<T, dim - 1>)) )), 
    buffer_size(1) {
    new ((char*)buffer) Grid<T, dim - 1>(t);
}

template <typename T>
Grid<T, 1>::Grid(T const &t) :
    //конструктор (неявное преобразование типов T -> Grid<T, 1>)

    Grid<T, 1>(1, t) { };

template <typename T, std::size_t dim>
Grid<T, dim>::Grid() : 
    //конструктор по умолчанию (создает пустую сетку)

    buffer(nullptr), 
    buffer_size(0) { }

template <typename T>
Grid<T, 1>::Grid() : 
    //конструктор по умолчанию (создает пустую сетку)

    buffer(nullptr), 
    buffer_size(0) { }

template <typename T, std::size_t dim>
Grid<T, dim>::Grid(Grid<T, dim> const &copy) :
    //конструктор копирования 

    buffer(reinterpret_cast<Grid<T, dim - 1>*>( operator new(sizeof(Grid<T, dim - 1>) * copy.get_dim_size(0)) )), 
    buffer_size(copy.get_dim_size(0)) {
    for (size_type idx = 0; idx != buffer_size; ++idx) {
        new ((char*)buffer + idx * sizeof(Grid<T, dim - 1>)) Grid<T, dim - 1>(copy[idx]);
    }
}

template <typename T>
Grid<T, 1>::Grid(Grid<T, 1> const &copy) :
    //конструктор копирования 
    
    buffer(reinterpret_cast<T*>( operator new(sizeof(T) * copy.get_dim_size(0)) )), 
    buffer_size(copy.get_dim_size(0)) {
    for (std::size_t idx = 0u; idx != buffer_size; ++idx) {
        new ((char*)buffer + idx * sizeof(T)) T(copy[idx]);
    }
}

template <typename T, std::size_t dim>
Grid<T, dim>::Grid(Grid<T, dim> &&copy) : 
    //конструктор перемещения

    buffer(copy.buffer), 
    buffer_size(copy.buffer_size) {
    copy.buffer = nullptr;
    copy.buffer_size = 0;
}

template <typename T>
Grid<T, 1>::Grid(Grid<T, 1> &&copy) : 
    //конструктор перемещения

    buffer(copy.buffer), 
    buffer_size(copy.buffer_size) {
    copy.buffer = nullptr;
    copy.buffer_size = 0;
}

template <typename T, std::size_t dim>
Grid<T, dim>& Grid<T, dim>::operator=(Grid<T, dim> const &src) { 
    //оператор присваивания копированием

    if (std::addressof(src) == this) return *this;
    Grid<T, dim> tmp(src); 
    std::swap(buffer, tmp.buffer);
    std::swap(buffer_size, tmp.buffer_size);
    return *this;
}

template <typename T>
Grid<T, 1>& Grid<T, 1>::operator=(Grid<T, 1> const &src) { 
    //оператор присваивания копированием

    if (std::addressof(src) == this) return *this;
    Grid<T, 1> tmp(src); 
    std::swap(buffer, tmp.buffer);
    buffer_size = tmp.buffer_size;
    return *this;
}

template <typename T, std::size_t dim>
Grid<T, dim>& Grid<T, dim>::operator=(Grid<T, dim> &&src) {
    //оператор присваивания перемещением

    if (this == std::addressof(src)) return *this;
    Grid<T, dim> tmp(std::move(src));
    std::swap(this->buffer, tmp.buffer);
    buffer_size = tmp.buffer_size;
    return *this;
}

template <typename T>
Grid<T, 1>& Grid<T, 1>::operator=(Grid<T, 1> &&src) {
    //оператор присваивания перемещением

    if (this == std::addressof(src)) return *this;
    Grid<T, 1> tmp(std::move(src));
    std::swap(this->buffer, tmp.buffer);
    buffer_size = tmp.buffer_size;
    return *this;
}

template <typename T, std::size_t dim>
Grid<T, dim>::~Grid() { 
    // деструктор
    
    for(size_type idx = 0; idx != buffer_size; ++idx) {
        buffer[idx].~Grid<T, dim - 1>();
    }
    delete[] buffer;
}

template <typename T>
Grid<T, 1>::~Grid() { 
    // деструктор
    
    for(size_type idx = 0; idx != buffer_size; ++idx) {
        buffer[idx].~T();
    }
    delete[] buffer;
}


template <typename T, std::size_t dim>
typename Grid<T, dim>::size_type Grid<T, dim>::get_dim_size(std::size_t dim_num) const {
    //метод возвращает размер сетки по данному направлению

    if (buffer == nullptr) {
        std::cout << "Empty grid\n";
        return 0;
    }
    if (dim_num == 0) {
        return buffer_size;
    }
    return buffer[0].get_dim_size(dim_num - 1);
}

template <typename T>
typename Grid<T, 1>::size_type Grid<T, 1>::get_dim_size(std::size_t dim_num) const {
    //метод возвращает размер сетки по последнему направлению

    if (buffer == nullptr) {
        std::cout << "Empty grid\n";
        return 0;
    }
    if (dim_num == 0) {
        return buffer_size;
    }
    std::cout << "Dimension number is greater than the dimension of the grid\n";
    return 0;
}


template <typename T, std::size_t dim>
Grid<T, dim - 1>& Grid<T, dim>::operator[](size_type idx) const { 
    //оператор индексирования
    return  buffer[idx]; 
}

template <typename T>
T& Grid<T, 1>::operator[] (size_type idx) const {
    //оператор индексирования

    if (buffer == nullptr) {
        std::cout << "Empty grid\n";
        exit(1);
    }
    if (idx >= buffer_size) {
        std::cout << "Invalid index\n";
        exit(1);
    }
    return buffer[idx];
}

template <typename T>
T& Grid<T, 1>::operator[] (size_type idx) {
    // оператор индексирования не константный

    if (buffer == nullptr) {
        std::cout << "Empty grid\n";
        exit(1);
    }
    if (idx >= buffer_size) {
        std::cout << "Invalid index\n";
        exit(1);
    }
    return buffer[idx];
}

template <typename T, std::size_t dim>
template <typename ...Dims>
T Grid<T, dim>::operator()(size_type idx, Dims... other) const {
    //другой оператор индексирования 

    if (buffer == nullptr) {
        std::cout << "Empty grid\n";
        exit(1);
    }
    if (idx >= buffer_size) {
        std::cout << "Invalid index\n";
        exit(1);
    }
    return buffer[idx](other...);
}

template <typename T>
T Grid<T, 1>::operator()(size_type idx) const {
    //другой оператор индексирования 

    if (buffer == nullptr) {
        std::cout << "Empty grid\n";
        exit(1);
    }
    if (idx >= buffer_size) {
        std::cout << "Invalid index\n";
        exit(1);
    }
    return buffer[idx];
}

template <typename T, std::size_t dim>
template <typename ...Dims>
T& Grid<T, dim>::operator()(size_type idx, Dims... other) {
    //другой оператор индексирования 

    if (buffer == nullptr) {
        std::cout << "Empty grid\n";
        exit(1);
    }
    if (idx >= buffer_size) {
        std::cout << "Invalid index\n";
        exit(1);
    }
    return buffer[idx](other...);
}

template <typename T>
T& Grid<T, 1>::operator()(size_type idx) {
    //другой оператор индексирования 

    if (buffer == nullptr) {
        std::cout << "Empty grid\n";
        exit(1);
    }
    if (idx >= buffer_size) {
        std::cout << "Invalid index\n";
        exit(1);
    }
    return buffer[idx];
}

template <typename T, std::size_t dim>
void Grid<T, dim>::print() const {
    //печать сетки

    if (buffer == nullptr) {
        std::cout << "Empty grid\n";
        return;
    }
    for (size_type idx = 0; idx != buffer_size; ++idx) {
        std::cout << "[";
        buffer[idx].print();
        std::cout << "]\n";
    }
}

template <typename T>
void Grid<T, 1>::print() const {
    //печать сетки

    if (buffer == nullptr) {
        std::cout << "Empty grid\n";
        return;
    }
    for (size_type idx = 0; idx != buffer_size; ++idx) {
        std::cout << " " << buffer[idx] << " ";
    }
}