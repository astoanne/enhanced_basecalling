import os
from pathlib import Path
from datetime import datetime
from matplotlib import pyplot as plt
import streamlit as st
from remora import io
import pandas as pd
import numpy as np
import pod5
from app.accuracy import calculate_accuracy_error_qscore_from_md_and_cigar, compute_accuracy, pretty_print_acc
from app.bam_utils import extract_md_tag, generate_fastq, get_cigar_and_md_tag
from app.bpd_gradients import analyze_signal_gradients, analyze_signal_gradients_no_plot
from app.bpd_ruptures import analyze_ruptures_breakpoints
from app.bpd_utils import filter_indices_by_mismatch_regions, filter_signal_indices, match_breakpoints, mismatch_regions_with_pairwisealigner, plot_breakpoints_with_labels
from app.cigar_utils import extract_cigar_strings
from app.core import batch_edit_bam
from app.utils import diff_str, find_floor_index



def main(test_data_root,read_id,input_path,original_aligned_path,second_aligned_path,output_folder,reference_filepath,edited_path,output_path,edited_filename):


    
    # Original unaligned read
    pod5_dr = pod5.DatasetReader(test_data_root)
    bam_fh = io.ReadIndexedBam(input_path)
    pod5_read = pod5_dr.get_read(read_id)
    bam_read = bam_fh.get_first_alignment(read_id)
    io_read_unaligned = io.Read.from_pod5_and_alignment(pod5_read, bam_read)
    st.subheader("Original Unaligned Read")
    st.write(f"Input Path: `{input_path}`")
    st.write(f"Basecalls Length: {io_read_unaligned.seq_len}")
    # st.write(f"Reference Mapping Length: {io_read_unaligned.ref_seq_len}")
    # st.write(f"Reference Location: {io_read_unaligned.ref_reg}")

    #Original aligned read
    bam_fh = io.ReadIndexedBam(original_aligned_path)
    bam_read = bam_fh.get_first_alignment(read_id)

    io_read_original = io.Read.from_pod5_and_alignment(pod5_read, bam_read)
    st.subheader("Original Aligned Read")
    st.write(f"Original Aligned Path: `{original_aligned_path}`")
    st.write(f"Basecalls Length: {io_read_original.seq_len}")
    st.write(f"Reference Mapping Length: {io_read_original.ref_seq_len}")
    st.write(f"Reference Location: {io_read_original.ref_reg}")

    st.subheader("Original Read Accuracy")
    temp_file = f"{output_folder}/tempfile_for_alignment.fq"
    original_acc_temp_sam_file = f"{test_data_root}/tempfile_for_alignment.sam"
    generate_fastq(io_read_original.seq, temp_file, read_id, quality_char="I")
    command = f'minimap2/minimap2 -ax lr:hq "{reference_filepath}" "{temp_file}" > "{original_acc_temp_sam_file}"'
    print(f"Running command: {command}") 
    os.system(command)
    acc_dict=compute_accuracy(original_acc_temp_sam_file)
    # acc_dict=calculate_accuracy_error_qscore_from_md_and_cigar(io_read_original.full_align['cigar'],extract_md_tag(io_read_original))
    pretty_print_acc(acc_dict)

    # Second method aligned read
    bam_fh = io.ReadIndexedBam(second_aligned_path)
    bam_read = bam_fh.get_first_alignment(read_id)

    io_read_second = io.Read.from_pod5_and_alignment(pod5_read, bam_read)
    st.subheader("Second method Aligned Read")
    st.write(f"Second method Aligned Path: `{second_aligned_path}`")
    st.write(f"Basecalls Length: {io_read_second.seq_len}")
    st.write(f"Reference Mapping Length: {io_read_second.ref_seq_len}")
    st.write(f"Reference Location: {io_read_second.ref_reg}")

    ## demonstrate
    ##plot original signal
    original_signal = io_read_unaligned.dacs
    start_default = 4400
    end_default = 4600

    # Streamlit UI
    st.title("Breakpoint Detection: Partial Signal for Demonstration")

    # Add a range slider for selecting start and end
    start, end = st.slider(
        "Select the range of indices",
        min_value=0,
        max_value=len(original_signal) - 1,
        value=(start_default, end_default),
        step=1,
    )

    # Display selected range
    # st.write(f"Selected range: Start = {start}, End = {end}")

    # Extract the selected segment
    selected_segment = original_signal[start:end]

    # Plot the selected segment
    st.subheader("Original Signal Segment")
    fig, ax = plt.subplots()
    ax.plot(range(start, end), selected_segment)
    ax.set_title("Original Signal")
    ax.set_xlabel("signal index")
    ax.set_ylabel("dacs")
    st.pyplot(fig)

    ###############################################
    st.subheader("Breakpoints Detected by Original Method") 
    ##breakpoint with original method (query_to_signal)
    algorithm_breakpoints=io_read_unaligned.query_to_signal
    start_base=find_floor_index(algorithm_breakpoints,start)
    end_base=find_floor_index(algorithm_breakpoints,end)
    algorithm_breakpoints=algorithm_breakpoints[start_base+1:end_base+1]
    
    fig, ax = plt.subplots(figsize=(10, 5))

    # Plot the selected segment of the signal
    ax.plot(range(start, end), original_signal[start:end], label="Original Signal", color="blue")

    # Add vertical lines for breakpoints
    for breakpoint in algorithm_breakpoints:
        if start <= breakpoint < end:  # Only include breakpoints within the selected range
            ax.axvline(x=breakpoint, color="red", linestyle="--", label=f"Breakpoint: {breakpoint}")

    # Customize the plot
    ax.set_title("Original Signal with Breakpoints")
    ax.set_xlabel("dacs")
    ax.set_ylabel("signal index")
    ax.grid(True)

    # Render the plot in Streamlit
    st.pyplot(fig)

    # display data
    df = pd.DataFrame([algorithm_breakpoints], columns=[f"Breakpoint {i+1}" for i in range(len(algorithm_breakpoints))])
    st.dataframe(df)

    #############################################
    st.subheader("Breakpoints Detected by Ruptures(Method 1)") 
    ruptures_breakpoints=analyze_ruptures_breakpoints(original_signal,start,end)

    fig, ax = plt.subplots(figsize=(10, 5))

    # Plot the selected segment of the signal
    ax.plot(range(start, end), original_signal[start:end], label="Original Signal", color="blue")

    # Add vertical lines for breakpoints
    for breakpoint in ruptures_breakpoints:
        if start <= breakpoint < end:  # Only include breakpoints within the selected range
            ax.axvline(x=breakpoint, color="red", linestyle="--", label=f"Breakpoint: {breakpoint}")

    # Customize the plot
    ax.set_title("Original Signal with Breakpoints")
    ax.set_xlabel("dacs")
    ax.set_ylabel("signal index")
    ax.grid(True)

    # Render the plot in Streamlit
    st.pyplot(fig)


    df = pd.DataFrame([ruptures_breakpoints], columns=[f"Breakpoint {i+1}" for i in range(len(ruptures_breakpoints))])
    st.dataframe(df)

    #############################################
    st.subheader("Breakpoints Detected by Peaks in Gradient(Method 2)") 
    signal_gradient_breakpoints=analyze_signal_gradients(original_signal,start,end)

    fig, ax = plt.subplots(figsize=(10, 5))

    # Plot the selected segment of the signal
    ax.plot(range(start, end), original_signal[start:end], label="Original Signal", color="blue")

    # Add vertical lines for breakpoints
    for breakpoint in signal_gradient_breakpoints:
        if start <= breakpoint < end:  # Only include breakpoints within the selected range
            ax.axvline(x=breakpoint, color="red", linestyle="--", label=f"Breakpoint: {breakpoint}")

    # Customize the plot
    ax.set_title("Original Signal with Breakpoints")
    ax.set_xlabel("dacs")
    ax.set_ylabel("signal index")
    ax.grid(True)

    # Render the plot in Streamlit
    st.pyplot(fig)


    df = pd.DataFrame([signal_gradient_breakpoints], columns=[f"Breakpoint {i+1}" for i in range(len(signal_gradient_breakpoints))])
    st.dataframe(df)

    #############################
    # Define the tolerance window size
    

    st.subheader("Mismatches in Breakpoints")
    tolerance = st.slider( 
        "Adjust Tolerance",  # Slider label
        min_value=0,         # Minimum value
        max_value=200,       # Maximum value
        value=5,             # Default value
        step=1               # Step size
    )
    st.text(f"Tolerance: A mismatch within {tolerance} is treated as a match.")

    
    # Match breakpoints for each method
    matched_ruptures_algorithm, matched_ruptures_method, fp_ruptures, fn_ruptures = match_breakpoints(algorithm_breakpoints, ruptures_breakpoints, tolerance)
    matched_peaks_algorithm, matched_peaks_method, fp_peaks, fn_peaks = match_breakpoints(algorithm_breakpoints, signal_gradient_breakpoints, tolerance)
    # Plot the results
    plot_breakpoints_with_labels(matched_ruptures_algorithm, matched_ruptures_method, fp_ruptures, fn_ruptures, 'Ruptures(Method1)')
    plot_breakpoints_with_labels(matched_peaks_algorithm, matched_peaks_method, fp_peaks, fn_peaks, 'Peaks in Gradient(Method 2)')

    # Print "false positives" and "false negatives" as lists
    st.markdown("### Mismatch Breakpoints")
    st.markdown("We will revisit these points and improve basecalling accuracy")
    # False positives and negatives for Ruptures
    st.markdown("#### Mismatch in Ruptures Breakpoints(Method 1) and Original Breakpoints")
    st.markdown(f"**Mismatch Ruptures:** `{fp_ruptures}`")
    st.markdown(f"**Mismatch Original:** `{fn_ruptures}`")

    # False positives and negatives for Peaks
    st.markdown("#### Mismatch in Peaks in Gradient Breakpoints(Method 2) and Original Breakpoints")
    st.markdown(f"**Mismatch Peaks in Gradient:** `{fp_peaks}`")
    st.markdown(f"**Mismatch Original:** `{fn_peaks}`")
    to_edit_signal_indices = fp_ruptures + fn_ruptures + fp_peaks + fn_peaks
    to_edit_signal_indices=sorted(set(to_edit_signal_indices))

    st.markdown("---")
    
    apply_filter = st.checkbox("Skip matched regions", value=False)
    mismatched_regions=mismatch_regions_with_pairwisealigner(
        io_read_original.seq, 
        io_read_original.ref_seq, 
        io_read_original.query_to_signal, 
        mode="global",
        treat_unaligned_ends_as_mismatch=True
    )
    # Conditionally filter the indices based on checkbox
    if apply_filter:
        to_edit_signal_indices = filter_indices_by_mismatch_regions(to_edit_signal_indices,mismatched_regions)
    
    filter_interval = st.slider(
        "Set Proximity Skip Interval",  # Slider label
        min_value=0,                   # Minimum value
        max_value=200,                 # Maximum value
        value=5,                      # Default value
        step=1                         # Step size
    )
    st.text(f"Proximity Skip Interval: Points closer than {filter_interval} signal index units to each other will not be edited.")
    st.markdown(f"**Points to Revisit:** `{filter_signal_indices(to_edit_signal_indices,filter_interval)}`")
    st.markdown(f"**Total Number of Points to Revisit:** `{len(filter_signal_indices(to_edit_signal_indices,filter_interval))}`")

    ##################################
    st.title("Breakpoint Detection: Full Sequence")
    ##################!!###################
    tolerance_full = st.slider( 
        "Adjust Tolerance",  # Slider label
        min_value=1,         # Minimum value
        max_value=200,       # Maximum value
        value=5,             # Default value
        step=1               # Step size
    )
    st.text(f"Tolerance: A mismatch within {tolerance} is treated as a match.")
    
   # Add a range slider for selecting start and end
    start=0
    end=len(io_read_unaligned.dacs)

    # Display selected range
    # st.write(f"Selected range: Start = {start}, End = {end}")

    # Extract the selected segment
    selected_segment = original_signal[start:end]

    # Plot the selected segment
    with st.expander("Original Signal Segment"):
        fig, ax = plt.subplots()
        ax.plot(range(start, end), selected_segment)
        ax.set_title("Original Signal")
        ax.set_xlabel("signal index")
        ax.set_ylabel("dacs")
        st.pyplot(fig)

    ###############################################
    with st.expander("Breakpoints Detection Process"): 
        st.subheader("Breakpoints Detected by Original Method")
        ##breakpoint with original method (query_to_signal)
        algorithm_breakpoints=io_read_unaligned.query_to_signal
        start_base=find_floor_index(algorithm_breakpoints,start)
        end_base=find_floor_index(algorithm_breakpoints,end)
        algorithm_breakpoints=algorithm_breakpoints[start_base+1:end_base+1]
        

        # display data
        df = pd.DataFrame([algorithm_breakpoints], columns=[f"Breakpoint {i+1}" for i in range(len(algorithm_breakpoints))])
        st.dataframe(df)

        #############################################
        st.subheader("Breakpoints Detected by Ruptures(Method 1)")
        ruptures_breakpoints=analyze_ruptures_breakpoints(original_signal,start,end)


        df = pd.DataFrame([ruptures_breakpoints], columns=[f"Breakpoint {i+1}" for i in range(len(ruptures_breakpoints))])
        st.dataframe(df)

        #############################################
        st.subheader("Breakpoints Detected by Peaks in Gradient(Method 2)")
        signal_gradient_breakpoints=analyze_signal_gradients_no_plot(original_signal,start,end)


        df = pd.DataFrame([signal_gradient_breakpoints], columns=[f"Breakpoint {i+1}" for i in range(len(signal_gradient_breakpoints))])
        st.dataframe(df)

        #############################

        st.subheader("Mismatches in Breakpoints")
            # Match breakpoints for each method
        matched_ruptures_algorithm, matched_ruptures_method, fp_ruptures, fn_ruptures = match_breakpoints(algorithm_breakpoints, ruptures_breakpoints, tolerance_full)
        matched_peaks_algorithm, matched_peaks_method, fp_peaks, fn_peaks = match_breakpoints(algorithm_breakpoints, signal_gradient_breakpoints, tolerance_full)
        

        # Print "false positives" and "false negatives" as lists
        st.markdown("### Mismatch Breakpoints")
        st.markdown("We will revisit these points and improve basecalling accuracy")
        # False positives and negatives for Ruptures
        st.markdown("#### Mismatch in Ruptures Breakpoints(Method 1) and Original Breakpoints")
        st.markdown(f"**Mismatch Ruptures:** `{fp_ruptures}`")
        st.markdown(f"**Mismatch Original:** `{fn_ruptures}`")

        # False positives and negatives for Peaks
        st.markdown("#### Mismatch in Peaks in Gradient Breakpoints(Method 2) and Original Breakpoints")
        st.markdown(f"**Mismatch Peaks in Gradient:** `{fp_peaks}`")
        st.markdown(f"**Mismatch Original:** `{fn_peaks}`")

    
    # to_edit_signal_indices = fp_peaks + fn_peaks
    to_edit_signal_indices = fp_peaks + fn_peaks+fp_ruptures + fn_ruptures
    # fp_ruptures + fn_ruptures
    # + fp_peaks + fn_peaks ##TODO
    to_edit_signal_indices=sorted(set(to_edit_signal_indices))##TODO

    # intervals=[ (935, 945), (2520, 2525), (4055, 4065), (4445, 4455), (5350, 5360), (5475, 5480), (5550, 5555), (5585, 5590), (6270, 6275), (6285, 6290), (6370, 6380), (6720, 6725), (6740, 6745), (8470, 8490), (8495, 8500), (8505, 8510), (8515, 8520), (8965, 8970), (10160, 10165), (10180, 10190), (10425, 10430), (10470, 10480), (10570, 10575), (10615, 10625), (10950, 10975), (11000, 11015), (11020, 11060), (11285, 11290), (11330, 11350), (12525, 12530), (12560, 12570), (12585, 12595), (12715, 12720), (12730, 12750), (12795, 12805), (14125, 14170), (14260, 14285), (14290, 14295), (15015, 15030), (15295, 15335), (15345, 15355), (15360, 15365), (15390, 15400), (15705, 15725), (16585, 16590), (16595, 16605), (16615, 16630), (16865, 16870), (17000, 17005), (17015, 17020), (18265, 18295), (18300, 18310), (18320, 18335), (18410, 18420), (19610, 19620), (20470, 20485), (20500, 20505), (21905, 21910), (22395, 22400), (22625, 22650), (22675, 22680), (22690, 22695), (22700, 22710), (22870, 22875), (24660, 24665), (24940, 24950), (24970, 24975), (24980, 24995), (25055, 25070), (25080, 25095), (25300, 25305), (27430, 27435), (27635, 27645), (27800, 27805), (27825, 27835), (27840, 27860), (27865, 27870), (28140, 28145), (28305, 28310), (28320, 28325), (28350, 28355), (28380, 28415), (28560, 28565), (28685, 28700), (28785, 28795), (29360, 29370), (30185, 30200), (30235, 30245), (30420, 30455), (30465, 30480), (31235, 31250), (33310, 33325), (33410, 33415), (33430, 33435), (33805, 33815), (34265, 34290), (34305, 34310), (34340, 34355), (34475, 34495), (34520, 34535), (34560, 34565), (34600, 34615), (34855, 34865), (34870, 34880), (34895, 34900), (34915, 34920), (34930, 34935), (34940, 34945), (35005, 35035), (35050, 35065), (35075, 35080), (35120, 35165), (35245, 35255), (35290, 35315), (35325, 35330), (35335, 35340), (35395, 35400), (35480, 35505), (35550, 35565), (35620, 35625), (35635, 35640), (35650, 35660), (35685, 35690), (35695, 35700), (35715, 35725), (35780, 35785), (35840, 35870), (35905, 35910), (35915, 35930), (35965, 35980), (36000, 36010), (36050, 36070), (36150, 36155), (36215, 36225), (36230, 36240), (36260, 36280), (36325, 36335), (36345, 36350), (36460, 36470), (36510, 36515), (36560, 36575), (36580, 36590), (36595, 36615), (36670, 36675), (36690, 36725), (36820, 36830), (36850, 36865), (37025, 37035), (37155, 37160), (37225, 37230), (37495, 37505), (37510, 37515), (37565, 37570), (37620, 37630), (37690, 37705), (38095, 38105), (38560, 38570), (40400, 40410), (41280, 41290), (41300, 41320), (41340, 41345), (41480, 41495), (41520, 41525), (43255, 43290), (43300, 43305), (43320, 43325), (43345, 43355), (43360, 43365), (43370, 43380), (43385, 43425), (43460, 43475), (43495, 43505), (43555, 43570), (43730, 43740), (43800, 43815), (45050, 45055), (45205, 45685)]
    # to_edit_signal_indices = [(start + end) // 2 for start, end in intervals]
    if "filter_interval_full" not in st.session_state:
        st.session_state.filter_interval_full = 20
    # Initialize session state for the slider
    if "radius" not in st.session_state:
        st.session_state.radius = 5  # Default value
    if "apply_filter_full" not in st.session_state:
        st.session_state.apply_filter_full = True  # Default value for checkbox
      

    apply_filter_full = st.checkbox("Skip matched regions", value=st.session_state.apply_filter_full)
    st.session_state.apply_filter_full = apply_filter_full

    # Apply conditional logic based on checkbox
    if apply_filter_full:
        to_edit_signal_indices = filter_indices_by_mismatch_regions(to_edit_signal_indices, mismatched_regions)
    else:
        # Set session state values when filter is not applied
        if st.session_state.filter_interval_full != 60:  # Ensure we don't reset unnecessarily
            st.session_state.filter_interval_full = 60
        if st.session_state.radius != 7:  # Ensure we don't reset unnecessarily
            st.session_state.radius = 7
    
    filter_interval_full = st.slider(
        "Set Proximity Skip Interval",
        min_value=1,
        max_value=200,
        value=st.session_state.filter_interval_full,  # Bind to session state
        step=1,
        key="filter_interval_full"  # Specify the session state key
    )
    
    
    st.text(f"Proximity Skip Interval: Points closer than {filter_interval_full} signal index units to each other will not be edited.")
    st.markdown(f"**Total Number of Points to Revisit:** `{len(to_edit_signal_indices)}`")

    to_edit_signal_indices=filter_signal_indices(to_edit_signal_indices,filter_interval_full)
#################!!######################
#############start editing#############
    st.subheader("Sequence Edit")
    
    # Slider bound to session state
    radius = st.slider(
        "Set Radius",                 # Slider label
        min_value=3,                  # Minimum value
        max_value=50,                 # Maximum value
        value=st.session_state.radius,  # Bind to session state
        step=1,                       # Step size
        key="radius"                  # Specify the session state key
    )


    st.text(f"Radius: Includes {radius} bases around the target signal index for local sequence extraction and alignment.") 
    with st.expander("Original Sequence"): 
        st.text(io_read_unaligned.seq)
    st.markdown(f"**Original Sequence Length:** `{len(io_read_unaligned.seq)}`")
    new_seq1=batch_edit_bam(
        input_path=input_path,
        edited_path=edited_path,
        read_id=read_id,
        seq1=io_read_original.seq,
        seq2=io_read_second.seq,
        query_to_signal1=io_read_original.query_to_signal,
        query_to_signal2=io_read_second.query_to_signal,
        signal_indices=to_edit_signal_indices,
        radius=radius 
    )
    with st.expander("Final Edited Sequence"):
        st.markdown(f"{new_seq1}")
    st.markdown(f"**Edited Sequence Length:** `{len(new_seq1)}`")

    ########Re-aligned############
    st.subheader("Accuracy")
    temp_file = f"{output_folder}/tempfile_for_alignment.fq"
    generate_fastq(new_seq1, temp_file, read_id, quality_char="I")
    sam_file = f"{output_folder}/{edited_filename}"
    # Step 1: Run the command using the variable
    command = f'minimap2/minimap2 -ax lr:hq "{reference_filepath}" "{temp_file}" > "{sam_file}"'
    # command = f'~/dorado/bin/dorado aligner "{reference_filepath}" "{edited_path}" --output-dir "{output_folder}"'
    # command = f'~/minimap2/minimap2 -ax lr:hq "{reference_filepath}" "{edited_path}" > "{output_folder}/{edited_filename}"'
    # command = f'~/minimap2/minimap2 -a "{reference_filepath}" "{edited_path}" > "{output_folder}/test.sam"'
    print(f"Running command: {command}") 
    os.system(command)
    # st.markdown(f"Done. Realigned Edited BAM will be stored in:\n`realigned/{edited_filename}`")

    
    ####################################
    # edited_cigar, edited_md=get_cigar_and_md_tag(output_path,read_id)
    acc_dict_ori=compute_accuracy(original_acc_temp_sam_file)
    # acc_dict_new=calculate_accuracy_error_qscore_from_md_and_cigar(edited_cigar, edited_md)
    acc_dict_new=compute_accuracy(sam_file)
    pretty_print_acc(acc_dict_ori,acc_dict_new)
    with st.expander("CIGAR String Comparison"):
        # st.markdown(f"**Original CIGAR String:** `{io_read_original.full_align['cigar']}`")
        st.markdown(f"**Original CIGAR String:** `{extract_cigar_strings(sam_file)[0]}`")
        st.markdown(f"**New CIGAR String:** `{extract_cigar_strings(original_acc_temp_sam_file)[0]}`")
    st.write(f"Edited file stored in: `{edited_path}`, and re-alignment information stored in: `{sam_file}`")













    
    
#(matched_ruptures_algorithm, matched_ruptures_method, fp_ruptures, fn_ruptures), (matched_peaks_algorithm, matched_peaks_method, fp_peaks, fn_peaks)

    
if __name__ == "__main__":
    # port = int(os.environ.get("PORT", 8501))

    # Run the Streamlit app
    # os.system(f"streamlit run app.py --server.port={port} --server.address=0.0.0.0")

    st.title("Improved Basecalling with Breakpoint Detection")
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")

    edited_filename = f"edited-{timestamp}.bam"
    print("edited_bam will be stored in: "+edited_filename)
    # st.markdown(f"### Edited BAM will be stored in:\n`{edited_filename}`")


    # input_path = test_data_root / "dorado-basecalled-result-aligned-mv-tables-fast-unaligned.bam"
    # original_aligned_path = test_data_root / "dorado-basecalled-result-aligned-mv-tables-fast.bam"
    test_data_root = Path(".") / "data"
    input_path = test_data_root / "dorado-basecalled-result-unaligned-mv-tables-hac.bam"
    original_aligned_path = test_data_root / "dorado-basecalled-result-aligned-mv-tables-hac.bam"

    second_aligned_path = test_data_root / "aligned_can_map_sup.bam"

    output_folder = test_data_root / "realigned"
    reference_filepath = test_data_root / "chr13.mmi"

    edited_path =  test_data_root / edited_filename
    output_path = output_folder / edited_filename

    read_id = "fbf9c81c-fdb2-4b41-85e1-0a2bd8b5a138"

    main(test_data_root,read_id,input_path,original_aligned_path,second_aligned_path,output_folder,reference_filepath,edited_path,output_path,edited_filename)
