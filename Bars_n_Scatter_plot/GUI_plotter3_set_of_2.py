%reset -f

import tkinter as tk
from tkinter import messagebox
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random

def convert_to_numpy():
    try:
        text = text_box.get("1.0", tk.END).strip()
        if not text:
            messagebox.showwarning("Warning", "No data pasted!")
            return
        
        # Convert tab-separated text into a DataFrame
        data = [row.split('\t') for row in text.split('\n') if row]
        df = pd.DataFrame(data)
        
        # Convert all values to numeric, coercing errors to NaN and dropping them
        df = df.apply(pd.to_numeric, errors='coerce').dropna()
        
        # Convert DataFrame to NumPy array
        global data2plot
        data2plot = df.to_numpy()
        
        messagebox.showinfo("Success", "Data stored in variable 'data2plot'")
        
        # Close the GUI window
        root.destroy()
        
        # Plot the data
        plot_data(data2plot)
    except Exception as e:
        messagebox.showerror("Error", f"Failed to convert to numpy array: {e}")

def plot_data(data):
    # Convert to DataFrame for easier handling
    df = pd.DataFrame(data)
    total_cols = df.shape[1]
    num_pairs = total_cols // 2  # number of complete pairs
    has_extra = (total_cols % 2 == 1)  # True if there's an extra column

    if has_extra:
        # Pair only the first 2*num_pairs columns
        group1 = df.iloc[:, :2*num_pairs:2]
        group2 = df.iloc[:, 1:2*num_pairs:2]
        extra = df.iloc[:, -1]  # last column is extra
    else:
        group1 = df.iloc[:, ::2]
        group2 = df.iloc[:, 1::2]
    
    # Compute means and standard errors for the paired groups
    mean_group1 = group1.mean()
    mean_group2 = group2.mean()
    std_group1 = group1.std() / np.sqrt(group1.count())
    std_group2 = group2.std() / np.sqrt(group2.count())
    
    if has_extra:
        extra_mean = extra.mean()
        extra_std = extra.std() / np.sqrt(extra.count())
    
    # Set bar widths and positions for paired columns
    group_width = 0.6        # Total width for a pair of bars
    bar_width = group_width / 2  # Width of individual bars
    dot_spread = 0.05        # Horizontal spread for scatter dots
    
    x = np.arange(num_pairs)         # one x position per pair
    x1 = x - bar_width / 2           # positions for group 1 (left bars)
    x2 = x + bar_width / 2           # positions for group 2 (right bars)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Define Excel-like colors for paired bars (use first two colors)
    excel_colors = ['#4472C4', '#ED7D31', '#A5A5A5', '#FFC000', '#70AD47', '#FF0000']
    
    # Plot paired bars
    bars1 = ax.bar(x1, mean_group1, width=bar_width, color=excel_colors[0])
    bars2 = ax.bar(x2, mean_group2, width=bar_width, color=excel_colors[1])
    
    # Plot scatter points for paired groups
    for i, bar in enumerate(bars1):
        for measurement in group1.iloc[:, i].dropna():
            ax.plot(
                bar.get_x() + bar.get_width()/2 + random.uniform(-dot_spread, dot_spread),
                measurement,
                marker='o', markersize=4, color='black', linestyle='None'
            )
    for i, bar in enumerate(bars2):
        for measurement in group2.iloc[:, i].dropna():
            ax.plot(
                bar.get_x() + bar.get_width()/2 + random.uniform(-dot_spread, dot_spread),
                measurement,
                marker='o', markersize=4, color='black', linestyle='None'
            )
    
    # Add error bars for paired groups
    ax.errorbar(x1, mean_group1, yerr=std_group1, fmt='none',
                ecolor='red', elinewidth=5, capsize=8, capthick=3, zorder=3)
    ax.errorbar(x2, mean_group2, yerr=std_group2, fmt='none',
                ecolor='red', elinewidth=5, capsize=8, capthick=3, zorder=3)
    
    # If there is an extra column, plot it on its own to the right
    if has_extra:
        x_extra = np.array([num_pairs])  # place extra bar after the paired groups
        bar_extra = ax.bar(x_extra, extra_mean, width=bar_width, color=excel_colors[0])
        for measurement in extra.dropna():
            ax.plot(
                bar_extra[0].get_x() + bar_extra[0].get_width()/2 + random.uniform(-dot_spread, dot_spread),
                measurement,
                marker='o', markersize=4, color='black', linestyle='None'
            )
        ax.errorbar(x_extra, extra_mean, yerr=extra_std, fmt='none',
                    ecolor='red', elinewidth=5, capsize=8, capthick=4, zorder=5)
    
    # Remove all text elements except enlarged Y-axis tick labels
    ax.set_xticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    
    # Increase Y-axis tick label size
    ax.tick_params(axis='y', labelsize=20)
    
    plt.show()

# Create main window
root = tk.Tk()
root.title("Excel to NumPy Converter")

# Input text box
text_box = tk.Text(root, height=10, width=50)
text_box.pack(pady=10)

# Convert button
convert_button = tk.Button(root, text="Convert to NumPy", command=convert_to_numpy)
convert_button.pack()

# Run the application
root.mainloop()
