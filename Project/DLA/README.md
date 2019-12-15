# Hướng dẫn sử dụng mã nguồn
## Yêu cầu môi trường
1. Hệ điều hành *Ubuntu 18.04LTS*: Chương trình đã chạy thành công trên hệ điều hành *Ubuntu 18.04LTS* nhưng có thể gặp lỗi khi chạy trên các hệ điều hành khác.
2. Cài đặt *GCC* và *OpenMPI*.
3. Cài đặt *Python3* và *pip3*.

## Cài đặt
1. Mở Terminal và thay đổi thư mục làm việc đến thư mục project.
2. Chạy lệnh *pip3 install -r requirements.txt* để cài đặt các module Numpy và Matplotlib cho Python. Nếu không thực hiện thao tác này, chương trình vẫn chạy bình thường, nhưng không thể tạo ảnh và video.

## Thiết lập thông số và chạy chương trình
1. Có thể thiết lập các thông số trong file *config.txt*. Mỗi tham số được đặt trên 1 dòng, chi tiết được chú thích ở cuối file *config.txt* theo thứ tự.
2. Chạy lệnh *./dla.sh* để thực hiện 3 bước: Biên dịch, chạy mô phỏng và visualize. Có thể chạy từng bước bằng cách chạy từng lệnh trong file *dla.sh*.
3. Lưu ý: Với lệnh *mpirun -H localhost:5 -np 5 dla.exe*, giá trị *5* trong *localhost:5* chỉ số process tối đa mà máy tính ở địa chỉ localhost có thể chạy, có thể áp dụng tương tự cho các máy khác; giá trị *5* trong *-np 5* chỉ tổng số process sẽ chạy.
4. Kết quả mô phỏng sẽ được lưu vào thư mục *./log/*, *./image/* và *./video/*.