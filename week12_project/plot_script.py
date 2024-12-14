import numpy as np
import scipy.optimize as opt
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from colorsys import hsv_to_rgb as hsv
import pathlib
import sys

font = {'weight': 'normal',
        'size': 14}
matplotlib.rc('font', **font)

# Функция ввывода загрузки
def progress_bar(current, total, bar_length=50):
    fraction = current / total
    arrow = '#' * int(fraction * bar_length)
    spaces = ' ' * (bar_length - len(arrow))
    sys.stdout.write(f'\r[{arrow}{spaces}] {int(fraction * 100)}%')
    sys.stdout.flush()

print("\nLoading data into Python")

# Итерируемся по директориям
amp_path = pathlib.Path("./in_out/output/ampermeters")
vol_path = pathlib.Path("./in_out/output/voltmeters")

amp_files = list(amp_path.iterdir())
vol_files = list(vol_path.iterdir())

print("\nPlotting...")

total_files = len(amp_files) + len(vol_files)
processed_files = 0
progress_bar(processed_files, total_files)

# Графики с амперметров
for file in amp_files:
    plt.figure(figsize=[10, 7], dpi=300)
    df = pd.read_csv(str(file))
    plt.plot(df["t"][2:-1], df["i"][2:-1], color="red", lw=1)
    plt.xlabel("$Time$, sec")
    plt.ylabel("$Current$, A")
    plt.title("Ampermeter " + str(file)[-5])
    plt.savefig(str(file)[:-26] + "plots\\" + str(file)[-26:-3] + "png")
    plt.close()

    # Обновляем прогресс-бар
    processed_files += 1
    progress_bar(processed_files, total_files)

# Графики с вольтметров
for file in vol_files:
    plt.figure(figsize=[10, 7], dpi=300)
    df = pd.read_csv(str(file))
    plt.plot(df["t"][2:-1], df["v"][2:-1], color="blue", lw=1)
    plt.xlabel("$Time$, sec")
    plt.ylabel("$Voltage$, V")
    plt.title("Voltmeter " + str(file)[-5])
    plt.savefig(str(file)[:-25] + "plots\\" + str(file)[-25:-3] + "png")
    plt.close()

    # Обновляем прогресс-бар
    processed_files += 1
    progress_bar(processed_files, total_files)

print("\n\nVisualization is finished!")
