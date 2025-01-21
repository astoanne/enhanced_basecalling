import streamlit as st

# Set the title of the app
st.title("Simple Streamlit App")

# Add a text input widget
name = st.text_input("What's your name?", "Guest")

# Display a button and respond when it's clicked
if st.button("Say Hello"):
    st.write(f"Hello, {name}! Welcome to the Streamlit app!")
