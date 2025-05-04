import matplotlib.pyplot as plt

# Dane odkształcenia (strain) i naprężenia (stress)
x_values =  [0.000000000000000, 9.660103178023142, 20.035778612839660, 29.695881790862785, 29.695881790862785, 40.429328067990042, 49.731675689788460, 60.107351124604975]

y_values = [0.000000000000000, 13.933655262925175, 30.710904590835177, 40.663511153932475, 40.663511153932475, 45.924175056897390, 47.345975126684294, 45.071093192648561]

# Tworzenie wykresu
plt.plot(x_values, y_values, marker='o', linestyle='-', label="Krzywa naprężenie-odkształcenie")

# Opis osi i tytuł
plt.xlabel("Odkształcenie [%]")
plt.ylabel("Naprężenie [MPa]")
plt.title("Wykres naprężenie-odkształcenie")
plt.legend()
plt.grid(True, linestyle='--', linewidth=0.5)

# Wyświetlenie wykresu
plt.show()
