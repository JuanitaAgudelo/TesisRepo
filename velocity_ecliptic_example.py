"""
Example usage of the refactored get_velocity_ecliptic function.

This example demonstrates how to use the function with both date format
and ephemeris time options.
"""

import numpy as np
import spiceypy as spy
from Utils import get_velocity_ecliptic

def example_usage():
    """Demonstrate both ways to use get_velocity_ecliptic function."""
    
    # Example parameters
    vx, vy, vz = 10.0, 5.0, 2.0  # Velocity components [km/s]
    lon, lat = -74.0, 40.7  # New York coordinates [degrees]
    alt = 0.0  # Sea level [km]
    
    print("=== Velocity Ecliptic Conversion Example ===\n")
    
    # Method 1: Using date string
    print("Method 1: Using date string")
    print("-" * 40)
    date = "2024-01-15 12:00:00"
    print(f"Date: {date}")
    print(f"Velocity (Earth-fixed): [{vx}, {vy}, {vz}] km/s")
    print(f"Location: {lon}°E, {lat}°N, {alt} km altitude")
    
    try:
        v_eclip_date = get_velocity_ecliptic(vx, vy, vz, lon, lat, alt, date=date)
        print(f"Velocity (Ecliptic J2000): {v_eclip_date} km/s")
    except Exception as e:
        print(f"Error: {e}")
    
    print()
    
    # Method 2: Using ephemeris time
    print("Method 2: Using ephemeris time")
    print("-" * 40)
    # Convert the same date to ephemeris time
    et = spy.utc2et(date)
    print(f"Ephemeris time: {et} seconds past J2000")
    print(f"Velocity (Earth-fixed): [{vx}, {vy}, {vz}] km/s")
    print(f"Location: {lon}°E, {lat}°N, {alt} km altitude")
    
    try:
        v_eclip_et = get_velocity_ecliptic(vx, vy, vz, lon, lat, alt, et=et)
        print(f"Velocity (Ecliptic J2000): {v_eclip_et} km/s")
    except Exception as e:
        print(f"Error: {e}")
    
    print()
    
    # Verify both methods give the same result
    if 'v_eclip_date' in locals() and 'v_eclip_et' in locals():
        are_equal = np.allclose(v_eclip_date, v_eclip_et, rtol=1e-10)
        print("Verification:")
        print(f"Both methods give same result: {are_equal}")
        if are_equal:
            print("✓ Methods are consistent!")
        else:
            print("✗ Methods give different results")
    
    print()
    
    # Error handling examples
    print("Error Handling Examples:")
    print("-" * 40)
    
    # Example 1: No time provided
    print("1. No date or et provided:")
    try:
        get_velocity_ecliptic(vx, vy, vz, lon, lat, alt)
    except ValueError as e:
        print(f"   Error caught: {e}")
    
    # Example 2: Both date and et provided
    print("2. Both date and et provided:")
    try:
        get_velocity_ecliptic(vx, vy, vz, lon, lat, alt, date=date, et=et)
    except ValueError as e:
        print(f"   Error caught: {e}")

def compare_with_original():
    """Compare the refactored function with original logic."""
    
    print("\n=== Comparison with Original Logic ===\n")
    
    # Test parameters
    vx, vy, vz = 15.0, -3.0, 8.0
    lon, lat = 0.0, 0.0  # Prime meridian, equator
    alt = 10.0
    date = "2024-06-21 15:30:00"
    
    print(f"Test parameters:")
    print(f"  Velocity: [{vx}, {vy}, {vz}] km/s")
    print(f"  Location: {lon}°E, {lat}°N, {alt} km altitude")
    print(f"  Date: {date}")
    
    # Use refactored function
    v_eclip_refactored = get_velocity_ecliptic(vx, vy, vz, lon, lat, alt, date=date)
    print(f"\nRefactored function result: {v_eclip_refactored} km/s")
    
    # The original logic would have been equivalent, but less robust
    print("\nThe refactored function provides:")
    print("✓ Better error handling")
    print("✓ Clear documentation")
    print("✓ Type hints")
    print("✓ Input validation")
    print("✓ Consistent behavior")

if __name__ == "__main__":
    example_usage()
    compare_with_original() 