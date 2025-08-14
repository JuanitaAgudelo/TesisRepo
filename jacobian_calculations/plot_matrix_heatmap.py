import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_matrix_heatmap(matrix, title="6x6 Matrix Heatmap", figsize=(8, 6), 
                       cmap='viridis', annot=True, fmt='.2e', scientific_notation=True):
    """
    Plot a 6x6 matrix as a heatmap
    
    Parameters:
    -----------
    matrix : numpy.ndarray
        6x6 matrix to plot
    title : str
        Title for the plot
    figsize : tuple
        Figure size (width, height)
    cmap : str
        Colormap for the heatmap
    annot : bool
        Whether to annotate cells with values
    fmt : str
        Format for annotations ('.2e' for scientific notation, '.2f' for decimal)
    scientific_notation : bool
        Whether to use scientific notation for annotations
    """
    # Create figure and axis
    plt.figure(figsize=figsize)
    
    # Prepare annotation format
    if scientific_notation and annot:
        # Create custom annotations in scientific notation
        annot_data = np.array([[f"{val:.2e}" for val in row] for row in matrix])
        
        # Create heatmap with custom annotations
        sns.heatmap(matrix, 
                    annot=annot_data, 
                    fmt='', 
                    cmap=cmap, 
                    center=0,
                    square=True,
                    linewidths=0.5,
                    cbar_kws={"shrink": .8})
    else:
        # Create heatmap with standard formatting
        sns.heatmap(matrix, 
                    annot=annot, 
                    fmt=fmt, 
                    cmap=cmap, 
                    center=0,
                    square=True,
                    linewidths=0.5,
                    cbar_kws={"shrink": .8})
    
    plt.title(title, fontsize=14, fontweight='bold')
    plt.xlabel('Column Index', fontsize=12)
    plt.ylabel('Row Index', fontsize=12)
    plt.tight_layout()
    plt.show()

# Example usage
if __name__ == "__main__":
    # Example 1: Random 6x6 matrix
    np.random.seed(42)  # For reproducible results
    random_matrix = np.random.randn(6, 6)
    
    print("Example 1: Random 6x6 Matrix")
    print(random_matrix)
    plot_matrix_heatmap(random_matrix, 
                       title="Random 6x6 Matrix Heatmap",
                       cmap='coolwarm')
    
    # Example 2: Symmetric matrix
    symmetric_matrix = np.random.randn(6, 6)
    symmetric_matrix = (symmetric_matrix + symmetric_matrix.T) / 2
    
    print("\nExample 2: Symmetric 6x6 Matrix")
    print(symmetric_matrix)
    plot_matrix_heatmap(symmetric_matrix, 
                       title="Symmetric 6x6 Matrix Heatmap",
                       cmap='plasma')
    
    # Example 3: Identity-like matrix
    identity_like = np.eye(6) + 0.1 * np.random.randn(6, 6)
    
    print("\nExample 3: Identity-like 6x6 Matrix")
    print(identity_like)
    plot_matrix_heatmap(identity_like, 
                       title="Identity-like 6x6 Matrix Heatmap",
                       cmap='RdBu_r')
    
    # Example 4: Large values matrix (demonstrating scientific notation)
    large_matrix = np.random.randn(6, 6) * 1e12  # Very large values
    
    print("\nExample 4: Large Values Matrix (Scientific Notation)")
    print(large_matrix)
    plot_matrix_heatmap(large_matrix, 
                       title="Large Values Matrix (Scientific Notation)",
                       cmap='seismic',
                       scientific_notation=True)
    
    # Example 5: Your Jacobian matrix (if available)
    # If you have a jacobian matrix from your previous code:
    # plot_matrix_heatmap(jacobian, 
    #                    title="Jacobian Matrix Heatmap",
    #                    cmap='seismic',
    #                    scientific_notation=True) 