import os
import time
from random import randint
from shutil import copy
from zipfile import ZipFile
from os.path import basename
from flask import Flask, request, flash, redirect, send_from_directory
from flask import render_template
from werkzeug.utils import secure_filename

from processing_data2input.processing_update import add2file, basic_csv2format, split_gl, create_glstring, \
    add_races_from_manual_insertion
from run3parts import run_all

from apscheduler.schedulers.background import BackgroundScheduler


app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = os.path.abspath("files_from_server")


def remove_old_files():
    """
    remove file in 'app.config['UPLOAD_FOLDER']' that exist more than an hour
    :return: nothing
    """
    now = time.time()

    folder1 = app.config['UPLOAD_FOLDER']
    folder2 = os.path.abspath("static")

    files1 = [os.path.join(folder1, filename) for filename in os.listdir(folder1)]
    files2 = [os.path.join(folder2, filename) for filename in os.listdir(folder2)]
    files1.extend(files2)
    print(files1)
    for filename in files1:
        # if file exist more than hour, remove it (.zip, .png, .pdf)
        if (filename.endswith('zip') or filename.endswith('png') or filename.endswith('pdf')) and (now - os.stat(filename).st_mtime) > 3600:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except OSError:
                print("Unable to remove file: %s" % filename)


# every 5 hours, remove old files of users (.zip, .png, .pdf)
# can not remove them in regular way, because they were sent to user
sched = BackgroundScheduler(daemon=True)
sched.add_job(remove_old_files, 'interval', minutes=3000)
sched.start()

# global variables
submit_file = True
save_gl_file = False
select_f = False
select_m = False
child_count = 1
races_manually = ''

# for show the user details of the input he insert in gl string format:
# Accepted: father:A*02:01+A*30:04^B*40:01(..), mother:A*03:01+A*21:04^B*32:01+B*51:22^C*03(..)
gl_accepted = ""


@app.route('/', methods=['GET', 'POST'])
# @app.route('/home', methods=['GET', 'POST'])
def home():
    global submit_file
    global save_gl_file
    global child_count
    global select_f
    global select_m
    global gl_accepted
    global races_manually
    error_message = False
    famcode = ""
    is_serology = False
    race_dict = {}  # in case of insertion races to families. (key = FANCODE, value = list of races)

    if request.method == 'POST':
        try:
            if request.form["submit_btn"] == "Results":
                select_f, select_m = False, False  # release the blockage for choosing a parent that already chose
                gl_accepted = ''
                races_manually_temp = races_manually  # put the data in another variable and clear (use the tmp after)
                races_manually = ''
                is_serology = True if request.form["gentic_or_serology"] == 'Serology' else False

                if submit_file:  # the user submit a file
                    # check if the post request has the file part
                    if len(request.files) == 0:
                        flash('No file part')
                        return redirect(request.url)
                    # get first file
                    f = request.files['the_file']
                    # if user does not select file, browser also
                    # submit an empty part without filename
                    if f.filename == '':
                        return render_template('home.html', error_message=error_message,
                                               select_f=select_f, select_m=select_m, send2user=False,
                                               gl_accepted=gl_accepted,  races_manually=races_manually)

                    if f:
                        # before receive new file, remove all the files
                        # in order that won't be error from data that stay in system
                        for filename in os.listdir(app.config['UPLOAD_FOLDER']):
                            os.remove(os.path.join(app.config['UPLOAD_FOLDER'], filename))

                        filename = secure_filename(f.filename)
                        path2f_old = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                        f.save(path2f_old)  # todo: why save path2f_old and not path2f?
                        path2f = os.path.join(app.config['UPLOAD_FOLDER'], "final.csv")
                        # check if user inserted races, and add to race dict
                        # add_races_from_file_to_dict(path2f_old, race_dict)  # remove because new version
                        open_ambiguity = basic_csv2format(path2f_old, path2f, race_dict)
                        os.remove(path2f_old)

                if not submit_file:  # the user submit data manually, so create a file
                    submit_file = True  # change the value to the next insertion
                    # check if user inserted races, and add to race dict
                    add_races_from_manual_insertion(races_manually_temp, race_dict)
                    # convert file to a valid format
                    # maybe not need here, if I create from the start as valid format?
                    open_ambiguity = basic_csv2format(app.config['UPLOAD_FOLDER'] + "/data2file.csv", app.config['UPLOAD_FOLDER'] + "/final.csv", race_dict)
                    # os.remove(app.config['UPLOAD_FOLDER'] + "/data2file.csv")
                    path2f = app.config['UPLOAD_FOLDER'] + "/final.csv"
                    child_count = 1
                    save_gl_file = True  # flag for saving the file "data2file.csv", and add to zip for user

                # res_100 is a flag for grimm to run with 1000 results (can be true only if run_again is true)
                res_100 = False

                # run the code
                res_path, errors_path, run_again = run_all(path2f, ['A', 'B', 'C', 'DRB1', 'DQB1'],
                                                           app.config['UPLOAD_FOLDER'], res_100, is_serology,
                                                           race_dict, open_ambiguity)
                if run_again:
                    res_100 = True
                    res_path, errors_path, run_again = run_all(path2f, ['A', 'B', 'C', 'DRB1', 'DQB1'],
                                                               app.config['UPLOAD_FOLDER'], res_100, is_serology,
                                                               race_dict, open_ambiguity)
                race_dict.clear()  # need?

                # create random number for differ files of each user
                # in order for files not to be overwritten in case of parallel use by different users
                rand_user = str(randint(1, 101))

                # create zip for all the user files
                with ZipFile(''.join([app.config['UPLOAD_FOLDER'], "/output_to_user" + rand_user + ".zip"]),
                             "w") as zipObj:
                    zipObj.write(res_path, basename(res_path))
                    zipObj.write(errors_path, basename(errors_path))
                    for filename in os.listdir(app.config['UPLOAD_FOLDER']):
                        if filename.endswith(".png") or filename.endswith(".pdf"):
                            zipObj.write(os.path.join(app.config['UPLOAD_FOLDER'], filename),
                                         basename(os.path.join(app.config['UPLOAD_FOLDER'], filename)))

                    if save_gl_file:
                        os.rename(app.config['UPLOAD_FOLDER'] + "/data2file.csv",
                                  app.config['UPLOAD_FOLDER'] + "/user data in gl format.csv")
                        zipObj.write((app.config['UPLOAD_FOLDER'] + "/user data in gl format.csv"),
                                     basename(app.config['UPLOAD_FOLDER'] + "/user data in gl format.csv"))

                        save_gl_file = False  # initialize to next running

                old_name_p = os.path.join(app.config['UPLOAD_FOLDER'], "1.png")
                new_name_p = os.path.join(app.config['UPLOAD_FOLDER'], rand_user + "_1.png")
                # add the random number to the name of file "1.png" (if exist).
                # this file stay for visualization in website
                if os.path.exists(old_name_p):
                    os.rename(old_name_p, new_name_p)

                # move the file "{{rand_user}}_1.png" to static (for show to user in pop up window)
                if os.path.exists(new_name_p):
                    copy(new_name_p, os.path.join(os.path.abspath("static"), rand_user + "_1.png"))
                    os.remove(new_name_p)

                os.remove(path2f)
                os.remove(res_path)
                os.remove(errors_path)
                if os.path.isfile(app.config['UPLOAD_FOLDER'] + "/data2file.csv"):
                    os.remove(app.config['UPLOAD_FOLDER'] + "/data2file.csv")
                if os.path.isfile(app.config['UPLOAD_FOLDER'] + "/user data in gl format.csv"):
                    os.remove(app.config['UPLOAD_FOLDER'] + "/user data in gl format.csv")
                if os.path.isfile(app.config['UPLOAD_FOLDER'] + "/final.csv"):  # todo: remove the file 'final.csv'?
                    os.remove(app.config['UPLOAD_FOLDER'] + "/final.csv")
                for filename in os.listdir(app.config['UPLOAD_FOLDER']):
                    if filename.endswith(".png") or filename.endswith(".pdf"):
                        os.remove(os.path.join(app.config['UPLOAD_FOLDER'], filename))

                select_f, select_m = False, False

                return render_template('home.html', error_message=error_message, select_f=select_f,
                                       select_m=select_m, send2user=True, rand_user=rand_user, gl_accepted=gl_accepted,
                                       races_manually=races_manually)

            if request.form["submit_btn"] == "Add family member":  # add members family
                submit_file = False
                famcode = request.form["member"]
                if famcode == "Father":
                    select_f = True
                if famcode == "Mother":
                    select_m = True
                if request.form["gl_or_split"] == "gl":  # gl string
                    gl_str = request.form["gl_string"]
                    gl_accepted = ''.join([gl_accepted, famcode, ':', gl_str, ';\n'])
                    A1_data, B1_data, C1_data, DRB1_data, DQB1_data, A2_data, B2_data, C2_data, DRB2_data, DQB2_data = \
                        split_gl(gl_str)
                else:  # split insertion
                    A1_data = request.form["A1"]  # todo: create a function for these requests
                    B1_data = request.form["B1"]
                    C1_data = request.form["C1"]
                    DRB1_data = request.form["DRB1"]
                    DQB1_data = request.form["DQB1"]
                    A2_data = request.form["A2"]
                    B2_data = request.form["B2"]
                    C2_data = request.form["C2"]
                    DRB2_data = request.form["DRB2"]
                    DQB2_data = request.form["DQB2"]

                    # gl to show for user
                    cur_gl = create_glstring(A1_data, A2_data, B1_data, B2_data, C1_data, C2_data, DRB1_data, DRB2_data, DQB1_data, DQB2_data)
                    gl_accepted = ''.join([gl_accepted, famcode, ':', cur_gl, ' ; '])

                if famcode == "Child":
                    famcode = child_count
                    child_count += 1

                add2file(app.config['UPLOAD_FOLDER'], famcode, A1_data, B1_data, C1_data, DRB1_data,
                         DQB1_data, A2_data, B2_data, C2_data, DRB2_data, DQB2_data, temp_f='data2file.csv')
            if request.form["submit_btn"] == "Clear family":  # clear the data about family
                if os.path.isfile(app.config['UPLOAD_FOLDER'] + "/data2file.csv"):
                    os.remove(app.config['UPLOAD_FOLDER'] + "/data2file.csv")
                if os.path.isfile(app.config['UPLOAD_FOLDER'] + "/final.csv"):
                    os.remove(app.config['UPLOAD_FOLDER'] + "/final.csv")
                child_count = 1
                gl_accepted, races_manually, select_f, select_m = '', '', False, False
                return render_template('home.html', error_message=error_message, select_f=select_f,
                                       select_m=select_m, send2user=False, gl_accepted=gl_accepted,
                                       races_manually=races_manually)

            if request.form["submit_btn"] == "add_race":  # add race
                cur_race = request.form["races"]
                races_manually = cur_race if races_manually == '' else races_manually + ';' + cur_race
                return render_template('home.html', error_message=error_message, select_f=select_f,
                                       select_m=select_m, send2user=False, gl_accepted=gl_accepted,
                                       races_manually=races_manually)
        except:
            error_message = True
            save_gl_file = False
            is_serology = False  # not need?
            races_manually = ''  # ?
            # remove the files in files_from_server
            # in error, remove all the files in order that won't be another error from data that stay in system
            for filename in os.listdir(app.config['UPLOAD_FOLDER']):
                os.remove(os.path.join(app.config['UPLOAD_FOLDER'], filename))
    return render_template('home.html', error_message=error_message, select_f=select_f,
                           select_m=select_m, send2user=False, gl_accepted=gl_accepted, races_manually=races_manually)


@app.route("/help")
def help():
    return render_template("help.html")


@app.route("/example")
def example():
    return render_template("example.html")


@app.route("/about")
def about():
    return render_template("about.html")


@app.route("/races")
def races():
    return render_template("races.html")


@app.route("/download_example1")
def download_example1():
    return send_from_directory(directory="static", filename="example_file1.csv")


@app.route("/download_example2")
def download_example2():
    return send_from_directory(directory="static", filename="example_file2.csv")


@app.route("/download_example1_with_races")
def download_example1_with_races():
    return send_from_directory(directory="static", filename="example_file1_with_races.csv")


@app.route("/download_example2_with_races")
def download_example2_with_races():
    return send_from_directory(directory="static", filename="example_file2_with_races.csv")


@app.route("/download_output")
@app.route("/download_output/<rand_user>")
def download_output(rand_user):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename="output_to_user" + rand_user + ".zip", as_attachment=True)


# clear the cache (to update the visualization in each running)
@app.after_request
def add_header(r):
    """
    Add headers to both force latest IE rendering engine or Chrome Frame,
    and also to cache the rendered page for 10 minutes.
    """
    r.headers["Cache-Control"] = "no-cache, no-store, must-revalidate"
    r.headers["Pragma"] = "no-cache"
    r.headers["Expires"] = "0"
    r.headers['Cache-Control'] = 'public, max-age=0'
    return r


if __name__ == "__main__":
    app.run(debug=True)
    # app.run(debug=True, host='0.0.0.0')
