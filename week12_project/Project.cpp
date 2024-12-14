#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <math.h>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <unordered_map>

# define PI 3.1415926535897932384626433832795

// Удаляет все элементы из вектора

template<typename T>
void delete_vector_elements(std::vector<T>& vec) {
    for (auto& elem : vec) {
        delete elem;
    }
    vec.clear();
}

// Матрица - вектор векторов, её печать (для проверки)

void print_matrix(std::vector<std::vector<double>> matrix) {
    for (int i = 0; i != matrix.size(); ++i) {
        for (int j = 0; j != matrix[i].size(); ++j) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

// Функция печати загрузки

void loading(int percent, int len = 50) {
    int position = (percent * len) / 100;

    std::cout << "\r[";

    for (int i = 0; i < len; ++i){
        if (i < position)
            std::cout << "#";
        else
            std::cout << " ";
    }

    if( percent < 10 ) 
        std::cout << "]   " << percent << "% ";
    else if ( percent < 100 )
        std::cout << "]  " << percent  << "% ";
    else
        std::cout << "] " << percent  << "% ";

    std::cout.flush();
}

// Решение системы алгиритмом Гаусса, in - матрица системы, out - вектор значений переменных

std::vector<double> gauss_solve(std::vector<std::vector<double>> matrix) {
    size_t checksum = 0;
    size_t matrix_height = matrix.size();
    size_t matrix_width = matrix_height + 1; // так как квадратная матрица + неоднородность СЛУ
    std::vector<double> answer;

    //Проверка, что входной вектор векторов - матрица

    for (size_t i = 0; i != matrix_height; ++i) {
        if (matrix[i].size() != matrix_width) { ++checksum; };
    }
    assert(checksum == 0 && "Wrong matrix dimensions!");

    //Применение алгоритма Гаусса

    //Прямое направление 

    for (size_t j = 0; j != matrix_height - 1; ++j) {

        double coef = matrix[j][j];

        // Ищем строку с j-тым коэффициентом не 0 и переставляем с исходной 

        std::vector<double> tmp; // для переставления
        size_t counter = j; // ищем

        while (coef == 0) {
            ++counter;
            assert(counter != matrix_height && "System has 0 or infinite number of solutuons"); // если не нашли
            coef = matrix[counter][j];
        }

        // переставили
        tmp = matrix[j];
        matrix[j] = matrix[counter];
        matrix[counter] = tmp;

        // Зануляем j-ые коэффициенты ниже

        for (size_t i = j + 1; i != matrix_height; ++i) {
            double mul = matrix[i][j] / coef;

            for (size_t k = j; k != matrix_width; ++k) {
                matrix[i][k] -= mul * matrix[j][k];
            }
        }
    }

    // Теперь обратное направление наверх

    for (size_t j = matrix_height - 1; j != 0; --j) {
        double coef = matrix[j][j];

        assert(coef != 0 && "System has 0 or infinite number of solutuons");

        for (size_t i = j - 1; i != -1; --i) {
            double mul = matrix[i][j] / coef;

            for (size_t k = j; k != matrix_width; ++k) {
                matrix[i][k] -= mul * matrix[j][k];
            }
        }
    }

    // Получаем ответ

    for (size_t j = 0; j != matrix_height; ++j) {
        answer.push_back(matrix[j][matrix_width - 1] / matrix[j][j]);
    }

    return answer;
}

// Node -- структура узлов, Bar -- абстрактных элементов (перемычек)
// Node и Bar делались с возможностью добавления элементам полярностей, то есть добавлением таких элементов как диоды, но в данной работе это не было реализованно

struct Node;

struct Bar {
private:
    double cur = 0; // ток 
    double vol = 0; // напряжение
    // между какими подключена
    Node* startNode;
    Node* endNode;
    std::string name; // имя

public:

    Bar(std::string name_) {
        startNode = nullptr;
        endNode = nullptr;
        name = name_;
    }

    Node* get_startNode() {
        return startNode;
    }

    Node* get_endNode() {
        return endNode;
    }

    Bar(Node* start, Node* end) {
        startNode = start;
        endNode = end;
    }

    void set_vol(double vol_) {
        vol = vol_;
    }

    void set_cur(double cur_) {
        cur = cur_;
    }

    double get_vol() {
        return vol;
    }

    double get_cur() {
        return cur;
    }

    // для вывода имен при проверки корректности ввода, будет в main
    virtual void print_name() {
        std::cout << name << std::endl;
    }

    // присоединение к узлу слева
    void connect_start(Node* node) {
        this->startNode = node;
    }

    // присоединение к узлу справа
    void connect_end(Node* node) {
        this->endNode = node;
    }

    // как определяются параметры элемента в зависимости от времени и итерации
    virtual void vac(double time, double step) {}


    virtual char get_type() {
        return 'B';
    }

    virtual double get_res() {
        return 0;
    }

    virtual char get_tag() { // отличает обычные элементы от осцилографов, у обычных элементов 'B', у осцилографов будут 'V' и 'A'
        return 'B';
    }

    virtual void output(std::ofstream*, double time) {
    }
};

// Узлы

struct Node {
private:
    std::vector<Bar*> ins; // перемычки, присоединенные к узлу endом
    std::vector<Bar*> outs; // перемычки, присоединенные к узлу startом
    unsigned int in_num = 0; // количество первых перемычек
    unsigned int out_num = 0; // и вторызх
    double potential; // потенциал узла, ключевой параметр

public:
    void set_potential(double potential_) {
        potential = potential_;
    }

    double get_potential() {
        return potential;
    }

    std::vector<Bar*> get_ins() {
        return ins;
    }

    std::vector<Bar*> get_outs() {
        return outs;
    }

    // присоединение элемента(премычки) к узлу с учётом направления
    void connect(Bar* element, char side) { // side - o(out) или i(in)
        if (side == 'o') {
            this->outs.push_back(element);
            element->connect_start(this); // элемент же подключается к узлу началом
            ++out_num;
        }
        else {
            this->ins.push_back(element);
            element->connect_end(this);
            ++in_num;
        }
    }

    // для проверки корректности ввода, будет в main
    void print() {
        std::cout << "This is a node with " << in_num << " ins and " << out_num << " outs.\n" << "ins:" << std::endl;
        for (size_t i = 0; i != ins.size(); ++i) {
            ins[i]->print_name();
        }
        std::cout << "outs:" << std::endl;
        for (size_t i = 0; i != outs.size(); ++i) {
            outs[i]->print_name();
        }
    }
};

// Резистор -- перемычка с сопротивлением, принципиально поменяется функция vac

struct Resistor : Bar
{
private:
    double res; // перемычка с сопротивлением

public:
    Resistor(double res_) : Bar("resistor") {
        res = res_;
    }

    Resistor(double res_, Node* start, Node* end) : Bar("resistor") {
        res = res_;
        start->connect(this, 'o'); // меняет и параметры Node и резистора
        end->connect(this, 'i');
    }

    void vac(double time, double step) { // для резистора параметры не будут зависить от времени, так как активный элемент
        this->set_vol(this->get_startNode()->get_potential() - this->get_endNode()->get_potential()); // напряжение как разница потенциалов
        this->set_cur(this->get_vol() / this->res); // для тока ещё поделим на сопротивление
    }

    double get_res() {
        return res;
    }

    char get_type() {
        return 'R';
    }

    void print_name() {
        std::cout << "resistor" << std::endl;
    }
};

// Источник постоянного напряжения -- перемычка с всегда заданным током через себя

struct Current_source : Bar
{
public:
    Current_source(double cur) : Bar("current source") {
        this->set_cur(cur);
    }

    Current_source(double cur, Node* start, Node* end) : Bar("current source") {
        this->set_cur(cur);
        start->connect(this, 'o');
        end->connect(this, 'i');
    }

    char get_type() {
        return 'I';
    }

    void vac(double time, double step) {
        this->set_vol(this->get_startNode()->get_potential() - this->get_endNode()->get_potential()); // напржение как разница потенциалов, ток всегда фиксированный
    }
};

// Источник напряжения как источник тока, но с фиксированным напряжением

struct Voltage_source : Bar
{
public:
    Voltage_source(double vol) : Bar("voltage source") {
        this->set_vol(vol);
    }

    Voltage_source(double vol, Node* start, Node* end) : Bar("voltage source") {
        this->set_vol(vol);
        start->connect(this, 'o');
        end->connect(this, 'i');
    }

    char get_type() {
        return 'V';
    }

    void vac(double time, double step) {}
};

// Синусоидальный источник тока -- тот у которого меняется параметр напряжение со временем

struct Sin_voltage_source : Voltage_source
{
private:
    double phase, freq, ampl;
public:
    // изначально создаётся в момент нулевой фазы и напряжения
    Sin_voltage_source(double ampl_, double freq_) : Voltage_source(0) {
        phase = 0;
        ampl = ampl_;
        freq = freq_;
    }

    // конструктор с сразу конектором
    Sin_voltage_source(double ampl_, double freq_, Node* start, Node* end) : Voltage_source((0), start, end) {
        phase = 0;
        ampl = ampl_;
        freq = freq_;
    }

    char get_type() {
        return 'V';
    }

    // тут уже возникает время, которое определет напряжение в конкретный момент
    void vac(double time, double step) {
        this->set_vol(ampl * sin(phase));
        phase = 2 * PI * time * freq;
    }
};

// Прямоугольник аналогично синусу только с другой функций vac от времени

struct Square_voltage_source : Voltage_source
{
private:
    double phase, freq, ampl;
public:
    Square_voltage_source(double ampl_, double freq_) : Voltage_source(0) {
        phase = 0;
        ampl = ampl_;
        freq = freq_;
    }

    Square_voltage_source(double ampl_, double freq_, Node* start, Node* end) : Voltage_source((0), start, end) {
        phase = 0;
        ampl = ampl_;
        freq = freq_;
    }

    char get_type() {
        return 'V';
    }

    // не синус, а прямоугольник
    void vac(double time, double step) {
        if (sin(phase) > 0) {
            this->set_vol(ampl);
        }
        else {
            this->set_vol(-ampl);
        }
        phase = 2 * PI * time * freq;
    }
};

// Провод -- источник напряжения с нулевым напряжением

struct Wire : Voltage_source
{
public:
    Wire() : Voltage_source(0) {
    }

    Wire(Node* start, Node* end) : Voltage_source(0, start, end) {
        this->set_vol(0);
    }

};

// Конденсатор -- источник переменного напряжения, для получения которого мы итерационно интегрируем

struct Capacitor : Voltage_source
{
private:
    double cap; // ёмкость
    double charge_0;
    double charge;  // напряжение будет определять заряд

public:
    Capacitor(double cap_, double vol_0_) : Voltage_source(0) { // начальные параметры задаём как ёмкость и начальное напряжение, т.к. оно более удобно для ввода
        cap = cap_;
        charge_0 = vol_0_ * cap;
        charge = charge_0;
    }

    // конструктор с коннектором
    Capacitor(double cap_, double vol_0_, Node* start, Node* end) : Voltage_source(0, start, end) {
        cap = cap_;
        charge_0 = vol_0_ * cap;
        charge = charge_0;
    }

    // тут уже влияют шаги по времени, т.к. итерационно интегрируем
    void vac(double time, double step) {
        charge -= this->get_cur() * step;  //  dq = - I * dt
        this->set_vol(charge / cap);
    }

    void print_name() {
        std::cout << "capacitor" << std::endl;
    }

};

// Индуктивность -- аналогично конденсатору, но источник переменного тока, который тоже получаем интегрированием

struct Inductor : Current_source
{
private:
    double ind; // индуктивность
    double cur_0; // начальный ток
public:
    Inductor(double ind_, double cur_0_) : Current_source(cur_0_) {
        ind = ind_;
        cur_0 = cur_0_;
        this->set_cur(cur_0);
    }

    // конструктер с коннектором
    Inductor(double ind_, double cur_0_, Node* start, Node* end) : Current_source(cur_0_, start, end) {
        ind = ind_;
        cur_0 = cur_0_;
        this->set_cur(cur_0);
    }

    char get_type() {
        return 'I';
    }

    void vac(double time, double step) {
        double cur_ = this->get_cur();
        this->set_vol(this->get_startNode()->get_potential() - this->get_endNode()->get_potential());  // напряжение как разница потенциалов
        this->set_cur(cur_ + this->get_vol() * step / ind); // дискретно интегрируем, чтобы получить ток, т.к. dI = U * dt
    }

    void print_name() {
        std::cout << "inductor" << std::endl;
    }

};


// Осцилографы: вольтметры и амперметры

// Вольтметр -> не идёт ток -> источник нулевого тока с функцией ввывода своих показаний в файл

struct Voltmeter : Current_source
{
public:
    Voltmeter() : Current_source(0) {}

    Voltmeter(Node* start, Node* end) : Current_source(0, start, end) {}

    // записываем данные в переданный файл в формате csv
    void output(std::ofstream *output_file, double time) {
        (*output_file) << time << "," << this->get_vol() << std::endl;
    }

    void print_name() {
        std::cout << "voltmeter" << std::endl;
    }

    double GetCurrent() {
        return this->get_cur();
    }

    char get_tag() { // отличает вольтметр от обычных элементов, у обычных элементов 'B'
            return 'V';
    }
};

// Идеальный Амперметр -> провод с функцией ввывода своих показаний в файл

struct Ampermeter : Wire
{
    public:
        Ampermeter() : Wire() {}

        Ampermeter(Node* start, Node* end) : Wire(start, end) {}

        // записываем данные в переданный файл в формате csv
        void output(std::ofstream *output_file, double time) {
            (*output_file) << time << "," << this->get_cur() << std::endl;
        }

        void print_name() {
            std::cout << "ampermeter" << std::endl;
        }

        double GetCurrent() {
            return this->get_cur();
        }

        char get_tag() { // отличает амперметр от обычных элементов, у обычных элементов 'B'
            return 'A';
        }
};

// Структура схемы, содержащей все элементы и узлы

struct Circuit {
private:
    // узлы
    std::vector<Node*> nodes;
    // обычные элементы
    std::vector<Bar*> resistors;
    std::vector<Bar*> cur_sources;
    std::vector<Bar*> wires; // провода и источники напряжений (так получилось, менять не хочу)
    // осцилографы
    std::vector<std::ofstream*> amp_outputs; // названия выходных файлов амперметров
    std::vector<std::ofstream*> vol_outputs; // названия выходных файлов вольтметров

public:

    std::vector<std::ofstream*> get_amp_outputs() {
        return amp_outputs;
    }

    std::vector<std::ofstream*> get_vol_outputs() {
        return vol_outputs;
    }

    // напечатать какие элементы есть в схеме, для проверки, будет в main
    void print() {
        std::cout << "This is a circuit with following elements" << std::endl;
        for (size_t i = 0; i != resistors.size(); ++i) {
            resistors[i]->print_name();
        }
        for (size_t i = 0; i != wires.size(); ++i) {
            wires[i]->print_name();
        }
        for (size_t i = 0; i != cur_sources.size(); ++i) {
            cur_sources[i]->print_name();
        }
        for (size_t i = 0; i != nodes.size(); ++i) {
            nodes[i]->print();
        }
    }

    // добавить узел
    void add_node(Node*& node) {
        nodes.push_back(node);
    }

    // добавить элемент (перемычку)
    void add_bar(Bar*& bar) {
        if (bar->get_type() == 'R') {
            resistors.push_back(bar);
        };

        if (bar->get_type() == 'V') {
            wires.push_back(bar);
        }

        if (bar->get_type() == 'I') {
            cur_sources.push_back(bar);
        }

        // Если осцилограф, добавляем к осцилографам и создаём название файла для сохранения данных с него

        // амперметр
        if (bar->get_tag() == 'A') {
            std::ostringstream oss;
            oss << "./in_out/output/ampermeters/output_A_" << amp_outputs.size() + 1 << ".csv"; // генерим название
            amp_outputs.push_back(new std::ofstream(oss.str()));
            oss.clear();
        }

        // вольтметр
        if (bar->get_tag() == 'V') {
            std::ostringstream oss;
            oss << "./in_out/output/voltmeters/output_V_" << vol_outputs.size() + 1 << ".csv";
            vol_outputs.push_back(new std::ofstream(oss.str()));
            oss.clear();
        }
    }

    // Одна из главных функций программы -- функция "решения" схемы. Так как подобный метод решения электрических схем это метод узловых потенциалов центральными являются они.
    
    // "Решить схему" значит найти все потенциалы узлов а так же ещё одни переменные и неизвестные параметры -- токи через источники напряжения (wires).

    // Тогда составляем систему уравнений относительно напряжений в узлах и токов через wires. В ней уравнения двух типов, первый тип - ура-ия правила Кирхгофа для всех* узлов, второй -- уравнения вида Voltage_wire_i = potential_Nodeend_i - potential_nodestart_i (определение источников токов)
    // * всех кроме первого, для перовго узла выбираем потенциал ноль. Это делается для того чтобы был потенциал для относительного отсчёта других и мы не получили линейно зависимую систему.

    void solve(double time, double step) {
        size_t height = nodes.size() + wires.size(); // количество неизвестных -- высота матрицы СЛУ, количество узлов + источников напряжения
        size_t width = nodes.size() + wires.size() + 1; // + неоднородность СЛУ -- ширины матрицы СЛУ

        std::vector<std::vector<double>> matrix;
        std::vector<double> tmp; // одна строка матрицы, всего их будет height строк

        // 1. Добавление уравнений первого типа в матрицу

        // 1.1. Установка потенциала первого узла на ноль во избежание неопределенности

        tmp.push_back(1);
        for (size_t i = 1; i != width; ++i) {
            tmp.push_back(0);
        }
        matrix.push_back(tmp);
        tmp.assign(width, 0);

        //  1.2 Добавление остальных nodes.size()-1 уравнений

        for (size_t i = 1; i != nodes.size(); ++i) {
            std::vector<Bar*> ins_ = nodes[i]->get_ins();
            std::vector<Bar*> outs_ = nodes[i]->get_outs();

            // 1.2.1 Элементы, указывающие внутрь рассматриваемого узла

            for (size_t j = 0; j < ins_.size(); ++j) {

                // 1.2.1.1 Токи от источников напряжения в уравнении

                if (ins_[j]->get_type() == 'V') {
                    size_t idx = nodes.size() + std::distance(wires.begin(), std::find(wires.begin(), wires.end(), ins_[j]));

                    tmp[idx] += 1;
                }

                // 1.2.1.2 Вклады от резисторов в уравнения вида: potential_Nodeend_i/R - potential_nodestart_i/R

                else if (ins_[j]->get_type() == 'R') {

                    size_t idx = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), ins_[j]->get_startNode()));

                    assert(idx != nodes.size() && "no such node in this circuit!");

                    tmp[idx] += 1 / ins_[j]->get_res();
                    tmp[i] -= 1 / ins_[j]->get_res();

                }

                // 1.2.1.3 Токи от источников токов в уравнении

                else if (ins_[j]->get_type() == 'I') {

                    // добавление к общему току этого узла, т.е. если есть источник тока то в правиле Кирхгофа в правой части равнество не нулю, а константе -I_source
                     
                    tmp[width - 1] -= ins_[j]->get_cur(); 
                }

            }

            // 1.2.2 Элементы, указывающие наружу узла (всё аналогично, только знаки мееняются)

            for (size_t j = 0; j < outs_.size(); ++j) {

                // 1.2.2.1 Токи от источников напряжения в уравнении

                if (outs_[j]->get_type() == 'V') {
                    size_t idx = nodes.size() + std::distance(wires.begin(), std::find(wires.begin(), wires.end(), outs_[j]));

                    tmp[idx] -= 1;
                }

                // 1.2.2.2 Вклады от резисторов в уравнения вида: potential_Nodeend_i/R - potential_nodestart_i/R

                else if (outs_[j]->get_type() == 'R') {

                    size_t idx = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), outs_[j]->get_endNode()));

                    assert(idx != nodes.size() && "no such node in this circuit!");

                    tmp[idx] += 1 / outs_[j]->get_res();
                    tmp[i] -= 1 / outs_[j]->get_res();
                }

                // 1.2.2.3 Токи от источников токов в уравнении
                   
                else if (outs_[j]->get_type() == 'I') {

                        // аналогично 1.2.1.3
                        tmp[width - 1] += outs_[j]->get_cur();
                    }
                }
             
            // 1.2.3 Итого сформировали строку, добавляем в матрицу

            matrix.push_back(tmp);
            tmp.assign(width, 0);
        }
        // к этому моменту есть node.size() уравнений

        // 2. Добавление уравнений второго типа в матрицу, тут всё проще

        // Проходимся по всем источникам напряжения
        for (size_t i = 0; i != wires.size(); ++i) {

            // Находим каким столбцам соответствуют клемы подключения
            size_t start_idx = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), wires[i]->get_startNode()));
            size_t end_idx = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), wires[i]->get_endNode()));

            // Записываем Voltage_wire_i = potential_Nodeend_i - potential_nodestart_i
            tmp[start_idx] -= 1;
            tmp[end_idx] += 1;
            tmp[width - 1] = wires[i]->get_vol();

            // Добавляем уравнение в систему
            matrix.push_back(tmp);
            tmp.assign(width, 0);
        }
        // итого к этому моменту имеем составленную матрицу СЛУ

        // Методом Гаусса решаем её!
        std::vector<double> solution = gauss_solve(matrix);

        // Присваиваем полученные потенциалы и токи узлам и wires соответственно
        for (size_t i = 0; i != nodes.size(); ++i) {
            nodes[i]->set_potential(solution[i]);
        }
        for (size_t i = 0; i != wires.size(); ++i) {
            wires[i]->set_cur(solution[i + nodes.size()]);
        }

        // Назначаем внутренние параметры элементов, зная все потенциалы и токи, это делает для каждого эоемента функция vac
        for (size_t i = 0; i != wires.size(); ++i) {
            wires[i]->vac(time, step);
        }
        for (size_t i = 0; i != resistors.size(); ++i) {
            resistors[i]->vac(time, step);
        }
        for (size_t i = 0; i != cur_sources.size(); ++i) {
            cur_sources[i]->vac(time, step);
        }
        // к этому моменту все параметры элементов (напряжение и ток через них и другие дополнительные) определены, осталось их только снять осцилографами

        // Снятие показаний схемы в данный момент осцилографами и их запись в выходные файлы

        // Вольтметры

        size_t osc_iterator = 0; // номер вольтметра

        for (size_t i = 0; i != cur_sources.size(); ++i) {
            if (cur_sources[i]->get_tag() == 'V') { // среди cur_sources нашли вольтметры
                cur_sources[i]->output(vol_outputs[osc_iterator], time); // запись в файл -- одна строка состояния схемы в данный момент
                ++osc_iterator;
            }
        }

        // Амперметры

        // аналогично
        osc_iterator = 0;

        for (size_t i = 0; i != wires.size(); ++i) {
            if (wires[i]->get_tag() == 'A') {
                wires[i]->output(amp_outputs[osc_iterator], time);
                ++osc_iterator;
            }
        }

        matrix.clear();
        tmp.clear();
    }

    virtual ~Circuit() {
        delete_vector_elements<Node*>(nodes);
        delete_vector_elements<Bar*>(resistors);
        delete_vector_elements<Bar*>(cur_sources);
        delete_vector_elements<Bar*>(wires);
        delete_vector_elements<std::ofstream*>(amp_outputs);
        delete_vector_elements<std::ofstream*>(vol_outputs);
    }
};

// Класс симуляции в целом, более удобная обёртка над цепью

struct Simulation
{
private:
    double total_time;
    double step;
    double time = 0;
    Circuit* circuit;
    std::vector<std::ofstream*> amp_outputs;
    std::vector<std::ofstream*> vol_outputs;

public:
    Simulation(Circuit* circuit_, double total_time_, double step_) {
        total_time = total_time_;
        step = step_;
        circuit = circuit_;
        amp_outputs = circuit->get_amp_outputs();
        vol_outputs = circuit->get_vol_outputs();
    }

    void run() {
        for (size_t i = 0; i < amp_outputs.size(); ++i) { // итерируемся по всем файлам для амперметром
            (*amp_outputs[i]) << "t,i" << std::endl; // вводим заголовоки в csv
        }

        for (size_t i = 0; i < vol_outputs.size(); ++i) { // аналогично
            (*vol_outputs[i]) << "t,v" << std::endl; 
        }

        double local_time = 0;
        size_t progress = 0;

        std::cout << "\nRunning simulation..." << std::endl;

        loading(progress);

        for (time = 0; time <= total_time; time += step) { // итерируемся по времени симуляции
            circuit->solve(time, step); 
            local_time += step;

            // печать прогресса симуляции в виде загрузки
            if (local_time > total_time * 0.01) {
                ++progress;
                loading(progress);
                local_time = 0;
            }
        }
        std::cout << std::endl;

        // закрытие файлов
        for (size_t i = 0; i < amp_outputs.size(); ++i) {
            (*amp_outputs[i]).close();
        }
        for (size_t i = 0; i < vol_outputs.size(); ++i) {
            (*vol_outputs[i]).close();
        }
        time = 0;
    }


    virtual ~Simulation() {
        delete this->circuit;
    }
};

// Функция конверции чисел с приставками  

double ValueInput(std::string value){
    double coef = 1;
    if(value[value.size()-1] == 'n') {
        coef = 1e-9;
    } else if (value[value.size()-1] == 'u') {
        coef = 1e-6;
    } else if (value[value.size()-1] == 'm') {
        coef = 1e-3;
    } else if (value[value.size()-1] == 'k') {
        coef = 1e3;
    } else if (value[value.size()-1] == 'M') {
        coef = 1e6;
    } else {
        return std::stod(value);
    }

    value.pop_back();

    double number = std::stod(value);

    return coef*number;
}

// Конверция текстового описания цепи в объект класса цепей, а также определение времён для дальнейшей симуляции

void ConvertInput(std::ifstream& file, double& simtime, double& simstep, Circuit* circuit) {
    
    std::string line;
    bool itsthestart = true;
    std::map<std::string, Node*> node_map; // стурктура с узлами, связанными с их названиями из файла

    // Автоматическая оценка времени -- посмотрим на характерное время процессов в контуре, это может быть RC, L/R, sqrt(LC) или частота генератора

    bool Rexist = false, Cexist = false, Lexist = false, Generator = false;
    double r, l, c, f_min, f_max;

    // Цикл счёта элементов

    while (getline(file, line)) {
        
        // пропускаем текст до -----
        if(itsthestart){
            if(line != "-----") {
                continue;
            } else { 
                itsthestart = false; 
                continue;
            }
        }
        
        // пропускаем пустые строки
        if (line == "") { 
            continue;
        }

        std::istringstream ss(line);

        std::string node_start;
        std::string node_end;
        std::string type;

        ss >> node_start;

        // если закончили перечисление элементов, выходим из цикла
        if (node_start == "-") {
                break;
            }

        ss >> node_end >> type;
        node_end.pop_back(); // убрать ":"
        
        // Добавляем в цепь новые узлы, которых ещё не было

        if (node_map.find(node_start) == node_map.end()) {
            Node* new_node = new Node();
            node_map[node_start] = new_node;
            circuit->add_node(new_node);
        }
        if (node_map.find(node_end) == node_map.end()) {
            Node* new_node = new Node();
            node_map[node_end] = new_node;
            circuit->add_node(new_node);
        }

        // Добавляем элементы в цепь, считывание параметров, а именно их количество, зависит от элемента

        if (type == "Wire" || type == "W") {
            Bar* new_wire = new Wire(node_map[node_start], node_map[node_end]);
            circuit->add_bar(new_wire);
        }
        else if (type == "Resistor" || type == "R") {
            std::string value;
            ss >> value;
            Bar* new_resistor = new Resistor(ValueInput(value), node_map[node_start], node_map[node_end]);
            circuit->add_bar(new_resistor);
            
            // для автоматической оценки времени
            Rexist = true; 
            r = ValueInput(value); 
        }
        else if (type == "Capacitor" || type == "C") {
            std::string value1, value2;
            ss >> value1 >> value2;
            Bar* new_inductor = new Capacitor(ValueInput(value1), ValueInput(value2), node_map[node_start], node_map[node_end]);
            circuit->add_bar(new_inductor);

            // для автоматической оценки времени
            Cexist = true;
            c = ValueInput(value1);
        }
        else if (type == "Inductor" || type == "L") {
            std::string value1, value2;
            ss >> value1 >> value2;
            Bar* new_inductor = new Inductor(ValueInput(value1), ValueInput(value2), node_map[node_start], node_map[node_end]);
            circuit->add_bar(new_inductor);

            // для автоматической оценки времени
            Lexist = true;
            l = ValueInput(value1);
        }
        else if (type == "Voltage_source" || type == "V_source") {
            std::string value;
            ss >> value;
            Bar* new_voltage_source = new Voltage_source(ValueInput(value), node_map[node_start], node_map[node_end]);
            circuit->add_bar(new_voltage_source);
        }
        else if (type == "Current_source" || type == "I") {
            std::string value;
            ss >> value;
            Bar* new_voltage_source = new Current_source(ValueInput(value), node_map[node_start], node_map[node_end]);
            circuit->add_bar(new_voltage_source);
        }
        else if (type == "Sin_voltage_source" || type == "Sin_V") {
            std::string value1, value2;
            ss >> value1 >> value2;
            Bar* new_voltage_source = new Sin_voltage_source(ValueInput(value1), ValueInput(value2), node_map[node_start], node_map[node_end]);
            circuit->add_bar(new_voltage_source);

            // для автоматической оценки времени
            if(!Generator){
                Generator = true; f_max = ValueInput(value2); f_min = ValueInput(value2);
            } else if (ValueInput(value2) > f_max) {
                f_max = ValueInput(value2);
            } else if (ValueInput(value2) < f_min ) {
                f_min = ValueInput(value2);
            }

        }
        else if (type == "Square_voltage_source" || type == "Sqr_V") {
            std::string value1, value2;
            ss >> value1 >> value2;
            Bar* new_voltage_source = new Square_voltage_source(ValueInput(value1), ValueInput(value2), node_map[node_start], node_map[node_end]);
            circuit->add_bar(new_voltage_source);

            // для автоматической оценки времени
            if(!Generator){
                Generator = true; f_max = ValueInput(value2); f_min = ValueInput(value2);
            } else if (ValueInput(value2) > f_max) {
                f_max = ValueInput(value2);
            } else if (ValueInput(value2) < f_min ) {
                f_min = ValueInput(value2);
            }

        }
        else if (type == "Ampermeter" || type == "A") {
            Bar* new_ampermeter = new Ampermeter(node_map[node_start], node_map[node_end]);
            circuit->add_bar(new_ampermeter);
        }
        else if (type == "Voltmeter" || type == "V") {
            Bar* new_voltmeter = new Voltmeter(node_map[node_start], node_map[node_end]);
            circuit->add_bar(new_voltmeter);
        }
    } // добавили все элементы в цепь

    // Определение или считывание параметров времени

    if(getline(file, line)){ // непустая строка
        std::istringstream ss(line);
        std::string value1, value2;
        ss >> value1 >> value2;
        simtime = ValueInput(value1);
        simstep = ValueInput(value2);
        if (simtime != 0 && simstep != 0){ // если не 0, то берём введенные пользователем
            return;
        }
    }

    // Наши оценки времён: RC, L/R, sqrt(LC) или 1/частота генератора

    std::vector<double> times;

    if(Rexist && Lexist){
        times.push_back(l/r);
    }
    if(Rexist && Cexist) {
        times.push_back(r*c);
    }
    if(Lexist && Cexist) {
        times.push_back(sqrt(l*c));
    }
    if(Generator) {
        times.push_back(1/f_max);
        times.push_back(1/f_min);
    }

    // Если совсем никаких характерных врёмен нет, возьмём условно случайное
    if(times.size() == 0) {
        simtime = 1;
        simstep = 1/10000;
        return;
    }

    // Если какие-то характерные времена нашлись, общее время моделяции оценивам сверху, а шаг снизу
    simtime = *(std::max_element(times.begin(), times.end()))*3; // три условных колебания системы
    simstep = *(std::min_element(times.begin(), times.end()))/1000;
}

int main() {
    std::cout << "\nStart of the program" << std::endl;
    std::ifstream file("./in_out/input/input.txt"); // файл входа

    double simtime_;
    double simstep_;
    Circuit* c1 = new Circuit;

    std::cout << "\nLoading data fom input" << std::endl;

    if (file.is_open()) {
        ConvertInput(file, simtime_, simstep_, c1); // конвертируем вход
        file.close();
    }

    // Если есть желание проверить корректность сцитывание, можно включить функцию вывода всех элементов и узлов
    // c1->print();

    // Создаём и запускаем симуляцию

    Simulation* sim1 = new Simulation(c1, simtime_, simstep_);
    sim1->run();
    // к этому моменту все csv файлы сделаны
    std::cout << "\nSimulation is finished!" << std::endl;

    // Визуализация файлов, построение графиков

    // Вызов скрипта на Python

    std::cout << "\nPython initialization..." << std::endl;

    std::string filename = "./plot_script.py";
    std::string command = "python ";
    command += filename;
    system(command.c_str());

    delete sim1;

    std::cout << "\nEnd of the program!\n" << std::endl;
    
    return 0;
}